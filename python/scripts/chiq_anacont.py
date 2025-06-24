#!/usr/bin/env python3

import os
import numpy as np

from chiq import __version__ as version
from chiq.pade import Pade
from chiq.spm import SpM
from chiq.bse_toml import load_params_from_toml
from chiq.mpi import COMM_WORLD as comm


def read_qw(chi_q_file):
    """Read q-points and Matsubara frequencies from chi_q_eigen.dat file.

    Parameters
    ----------
    chi_q_file : str
        Path to chi_q_eigen.dat file containing eigenvalues of susceptibility.

    Returns
    -------
    qs : list of str
        List of q-points in format "XX.YY.ZZ".
    omegas : list of int
        List of Matsubara frequency indices.
    nelems : int
        Number of eigenvalues per (q,omega) point.
    """
    omegas = []
    qs = []
    first_omega = None
    first_q = None
    nelems = None
    with open(chi_q_file, "r") as file:
        for line in file:
            line = line.strip()
            if line.startswith("#"):
                continue
            words = line.split()
            if len(words) < 2:
                continue
            omega = int(words[0])
            q = words[1]
            if first_omega is None:  # first record
                first_omega = omega
                first_q = q
                nelems = (len(words) - 2) // 2

            if omega == first_omega:
                qs.append(q)
            if q == first_q:
                omegas.append(omega)
    return qs, omegas, nelems


def read_chi_qw(chi_q_file):
    """Read chi(q,omega) from chi_q_eigen.dat file.

    Parameters
    ----------
    chi_q_file : str
        Path to chi_q_eigen.dat file containing eigenvalues of susceptibility.

    Returns
    -------
    chi_qw : dict of numpy.ndarray
        Dictionary of chi(q,omega) for each q-point.
        Keys are q-points in format "XX.YY.ZZ".
        Values are numpy.ndarray of shape (nelems, nomega) containing chi(q,omega) for each element.
    omegas : list of int
        List of Matsubara frequency indices.
    T : float
        Temperature.
    """
    qs, omegas, nelems = read_qw(chi_q_file)
    nomega = len(omegas)
    chi_qw = {q: np.zeros((nelems, nomega), dtype=np.complex128) for q in qs}

    omega2index = {omega: i for i, omega in enumerate(omegas)}

    T = None
    with open(chi_q_file, "r") as file:
        for line in file:
            line = line.strip()
            if line.startswith("#Temperature:"):
                T = float(line.split(":")[1].strip())
                continue
            if line.startswith("#"):
                continue
            if line.strip() == "":
                continue

            words = line.split()
            omega = int(words[0])
            q = words[1]
            cqw = np.array(
                [
                    complex(float(words[2 * i + 2]), float(words[2 * i + 3]))
                    for i in range(nelems)
                ],
                dtype=np.complex128,
            )
            chi_qw[q][:, omega2index[omega]] = cqw

    return chi_qw, omegas, T


def main():
    import argparse

    mpisize = comm.Get_size()
    mpirank = comm.Get_rank()

    parser = argparse.ArgumentParser(
        description="Analytic continuation of chi(q, iwn) to chi(q,w).",
        add_help=True,
    )

    parser.add_argument("toml", type=str, help="TOML configuration file")
    _version_message = f"ChiQ version {version}"
    parser.add_argument("--version", action="version", version=_version_message)

    args = parser.parse_args()

    params = load_params_from_toml(args.toml)

    dict_post = params["post"]
    post_output_dir = dict_post["output_dir"]

    dict_anacont = params["anacont"]
    input_file = os.path.join(post_output_dir, dict_anacont["input_file"])
    output_dir = dict_anacont["output_dir"]
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    output_prefix = dict_anacont["output_prefix"]
    wmax = dict_anacont["wmax"]
    wmin = dict_anacont["wmin"]
    wnum = dict_anacont["wnum"]
    method = dict_anacont["method"]
    eta = dict_anacont["eta"]
    loglambda = dict_anacont["loglambda"]
    maxiter = dict_anacont["maxiter"]
    initial_mu = dict_anacont["initial_mu"]

    ws = np.linspace(wmin, wmax, wnum)

    if mpisize > 1:
        if mpirank == 0:
            chi_qw, omegas, T = read_chi_qw(input_file)
        else:
            chi_qw, omegas, T = None, None, None

        chi_qw = comm.bcast(chi_qw, root=0)
        omegas = comm.bcast(omegas, root=0)
        T = comm.bcast(T, root=0)
    else:
        chi_qw, omegas, T = read_chi_qw(input_file)

    iwns = (2 * np.pi * T * 1j) * np.array(omegas, dtype=np.complex128)
    qs = list(chi_qw.keys())
    qs.sort()

    local_qs = qs[mpirank::mpisize]

    for q in local_qs:
        chi_iw = chi_qw[q]

        if dict_anacont["print_chi_q_iw"]:
            out_file = os.path.join(output_dir, f"{dict_anacont['output_prefix_chi_q_iw']}_{q}.dat")
            with open(out_file, "w") as file:
                for iomega, omega in enumerate(omegas):
                    file.write(f"{omega}")
                    for ielem in range(chi_iw.shape[0]):
                        file.write(
                            f" {np.real(chi_iw[ielem, iomega])} {np.imag(chi_iw[ielem, iomega])}"
                        )
                    file.write("\n")
            print(f"{out_file} is saved")

        nelem = chi_iw.shape[0]
        chi_w = np.zeros((nelem, len(ws)), dtype=np.complex128)
        if method == "pade":
            for ielem in range(nelem):
                pade = Pade(iwns, chi_iw[ielem, :])
                chi_w[ielem, :] = pade.evaluate(ws + eta * 1j)
        else:
            for ielem in range(nelem):
                spm = SpM(
                    beta=1.0 / T,
                    wmax=wmax * dict_anacont["wmax_factor"],
                    matsubara_points=2 * np.array(omegas),
                    chi_iwn=chi_iw[ielem, :],
                    l1_coeff=10**loglambda,
                    max_iter=maxiter,
                    initial_mu=initial_mu,
                )
                chi_w[ielem, :] = spm.evaluate(ws)

        out_file = os.path.join(output_dir, f"{output_prefix}_{q}.dat")
        with open(out_file, "w") as file:
            for iomega, omega in enumerate(ws):
                file.write(f"{omega}")
                for ielem in range(nelem):
                    file.write(f" {np.real(chi_w[ielem, iomega])}")
                    file.write(f" {np.imag(chi_w[ielem, iomega])}")
                file.write("\n")
        print(f"{out_file} is saved")


if __name__ == "__main__":
    main()
