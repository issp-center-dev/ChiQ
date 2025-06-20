#!/usr/bin/env python3

import numpy as np

from chiq import __version__ as version
from chiq.pade import Pade


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

    parser = argparse.ArgumentParser(
        description="Analytic continuation of chi.",
        add_help=True,
    )

    parser.add_argument("chi_q", type=str, help="Chi(w,q) file")
    _version_message = f"ChiQ version {version}"
    parser.add_argument("--version", action="version", version=_version_message)
    parser.add_argument("--wmax", type=float, default=10.0, help="Maximum frequency")
    parser.add_argument("--wmin", type=float, default=0.0, help="Minimum frequency")
    parser.add_argument("--wnum", type=int, default=101, help="Number of frequencies")
    parser.add_argument("--eta", type=float, default=1e-5, help="Imaginary shift")

    args = parser.parse_args()

    ws = np.linspace(args.wmin, args.wmax, args.wnum)
    eta = args.eta

    chi_qw, omegas, T = read_chi_qw(args.chi_q)
    iwns = (np.pi * T * 1j) * np.array(omegas, dtype=np.complex128)
    for q in chi_qw:
        chi_iw = chi_qw[q]
        nelem = chi_iw.shape[0]
        chi_w = np.zeros((nelem, len(ws)), dtype=np.complex128)
        for ielem in range(nelem):
            pade = Pade(iwns, chi_iw[ielem, :])
            chi_w[ielem, :] = pade.evaluate(ws + eta * 1j)

        with open(f"chiq_w_{q}.dat", "w") as file:
            for iomega, omega in enumerate(ws):
                file.write(f"{omega}")
                for ielem in range(nelem):
                    file.write(f" {np.real(chi_w[ielem, iomega])}")
                    file.write(f" {np.imag(chi_w[ielem, iomega])}")
                file.write("\n")


if __name__ == "__main__":
    main()
