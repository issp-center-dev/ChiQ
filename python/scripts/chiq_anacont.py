#!/usr/bin/env python3

import os
import sys
import logging

import numpy as np

from chiq import __version__ as version
from chiq.pade import Pade


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

    parser = argparse.ArgumentParser(
        description='Analytic continuation of chi.',
        add_help=True,
    )

    parser.add_argument('chi_q', type=str, help='Chi(w,q) file')
    _version_message = f'ChiQ version {version}'
    parser.add_argument('--version', action='version', version=_version_message)
    parser.add_argument('--wmax', type=float, default=10.0, help='Maximum frequency')
    parser.add_argument('--wmin', type=float, default=0.0, help='Minimum frequency')
    parser.add_argument('--wnum', type=int, default=101, help='Number of frequencies')
    parser.add_argument('--eta', type=float, default=1e-5, help='Imaginary shift')

    args = parser.parse_args()

    ws = np.linspace(args.wmin, args.wmax, args.wnum)
    eta = args.eta

    chi_qw, omegas, T = read_chi_qw(args.chi_q)
    iwns = (np.pi * T * 1j) * np.array(omegas, dtype=np.complex128)
    for q in chi_qw:
        chi_iw = chi_qw[q]
        # with open(f"chiq_iw_{q}.dat", "w") as file:
        #     for iomega, omega in enumerate(iwns):
        #         file.write(f"{np.imag(omega)}")
        #         file.write(f" {np.real(chi_iw[0, iomega])}")
        #         file.write(f" {np.imag(chi_iw[0, iomega])}")
        #         file.write("\n")

        pade = Pade(iwns, chi_iw[0, :])
        chi_w = pade.evaluate(ws + eta * 1j)
        with open(f"chiq_w_{q}.dat", "w") as file:
            for iomega, omega in enumerate(ws):
                file.write(f"{omega}")
                file.write(f" {np.real(chi_w[iomega])}")
                file.write(f" {np.imag(chi_w[iomega])}")
                file.write("\n")


if __name__ == "__main__":
    main()
