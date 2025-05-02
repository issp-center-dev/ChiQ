#!/usr/bin/env python3

import os
import sys
import logging

import numpy as np
from more_itertools import divide

from chiq import __version__ as version
from chiq.pade import Pade
# from chiq import bse_toml

# from chiq.mpi import comm


def read_qw(chi_q_file):
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
            if first_omega is None: # first record
                first_omega = omega
                first_q = q
                nelems = len(words) - 2

            if omega == first_omega:
                qs.append(q)
            if q == first_q:
                omegas.append(omega)
    return qs, omegas, nelems


def read_chi_qw(chi_q_file):
    qs, omegas, nelems = read_qw(chi_q_file)
    nomega = len(omegas)
    chi_qw = {q: np.zeros((nelems, nomega), dtype=np.float64) for q in qs}

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
            chi_qw[q][:, omega2index[omega]] = np.array(words[2:], dtype=np.float64)

    return chi_qw, omegas, T



def main():
    import argparse

    # rank = comm.Get_rank()
    # size = comm.Get_size()

    parser = argparse.ArgumentParser(
        description='Analytic continuation of chi.',
        add_help=True,
    )

    parser.add_argument('chi_q', type=str, help='Chi(w,q) file')
    _version_message = f'ChiQ version {version}'
    parser.add_argument('--version', action='version', version=_version_message)

    args = parser.parse_args()

    # ws = np.linspace(0, 10, 101)
    ws = np.array([0.5])
    eta = 1e-5

    # logger = logging.getLogger("bse")
    # fmt = "%(asctime)s %(levelname)s %(name)s :%(message)s"
    # # logging.basicConfig(level=logging.DEBUG, format=fmt)
    # logging.basicConfig(level=logging.INFO, filename="chiq_anacont_{}.log".format(rank), format=fmt, filemode='w')

    chi_qw, omegas, T = read_chi_qw(args.chi_q)
    iwns = (np.pi * T * 1j) * np.array(omegas, dtype=np.complex128)
    iwns_fit = iwns
    # iwns_fit = np.linspace(0, iwns[-1], 101)
    for q in chi_qw:
        chi_iw = chi_qw[q]
        pade = Pade(iwns, chi_iw[0, :])
        chi_iw_fit = pade.evaluate(iwns_fit)
        chi_w = pade.evaluate(ws + eta * 1j)

        with open(f"chiq_iw_{q}.dat", "w") as file:
            for iomega, omega in enumerate(iwns):
                file.write(f"{np.imag(omega)}")
                file.write(f" {np.real(chi_iw[0, iomega])}")
                file.write(f" {np.imag(chi_iw[0, iomega])}")
                file.write("\n")

        with open(f"chiq_iw_fit_{q}.dat", "w") as file:
            for iomega, omega in enumerate(iwns_fit):
                file.write(f"{np.imag(omega)}")
                file.write(f" {np.real(chi_iw_fit[iomega])}")
                file.write(f" {np.imag(chi_iw_fit[iomega])}")
                file.write("\n")

        with open(f"chiq_w_{q}.dat", "w") as file:
            for iomega, omega in enumerate(ws):
                file.write(f"{omega}")
                file.write(f" {np.real(chi_w[iomega])}")
                file.write(f" {np.imag(chi_w[iomega])}")
                file.write("\n")



if __name__ == "__main__":
    main()
