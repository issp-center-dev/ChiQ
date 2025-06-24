#!/usr/bin/env python3

import numpy as np
import os
import argparse

from chiq.chiqw_eigen_path import ChiQWEigenPath

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as colors


def main():
    P = argparse.ArgumentParser()
    P.add_argument("file_qpath", help="path of q-path file")
    P.add_argument("dir_chiqw", nargs="?", default=".", help="directory of input files")
    P.add_argument("--input_file_prefix", default="chi_q_w", help="prefix of input file (e.g., 'chi_q_w' for 'chi_q_w__XX.YY.ZZ.dat')")
    P.add_argument("--ielem", type=int, default=0, help="index of element to plot")
    P.add_argument("--wmin", type=float, default=None, help="minimum frequency to plot")
    P.add_argument("--wmax", type=float, default=None, help="maximum frequency to plot")
    P.add_argument("--part", type=str, default="real", choices=["real", "imag"], help="part of susceptibility to plot")
    P.add_argument("--output_file_prefix", default="chi_q_w", help="base name of output file (e.g., 'chi_q_w' for 'chi_q_w.pdf')")
    P.add_argument(
        "--format",
        default="pdf",
        type=lambda s: s.split(","),
        help="""output format. multiple formats can be specified by comma-separated string (e.g., 'pdf,png').
        Formats supported by savefig() of matplotlib.pyplot are available.""",
    )
    P.add_argument(
        "--data_out", default=None,
        help="""set filename to write chi(q,omega) along the q-path as text file (e.g., 'chiq_w.dat').
        If not specified, data will not be written."""
    )

    args = P.parse_args()

    print("\nRunning", os.path.basename(__file__))
    print(args)

    # --------------------------------------------------------------------------
    # get data to plot
    E = ChiQWEigenPath(args.file_qpath)
    xarray = E.get_x()
    ws, chiqw_re, chiqw_im = E.get_chiqw_on_path(dir_chiqw=args.dir_chiqw,
                                                 file_prefix=args.input_file_prefix,
                                                 ielem=args.ielem,
                                                 wmin=args.wmin,
                                                 wmax=args.wmax)

    n_q, n_w = chiqw_re.shape
    print("# of q points", n_q)
    print("# of w points", n_w)

    # --------------------------------------------------------------------------
    # plot

    xticks, xlabels = E.get_xticks(latex=True)

    fig, ax = plt.subplots(figsize=(8, 6))

    # Create meshgrid for pcolormesh
    X, Y = np.meshgrid(xarray, ws)

    # Plot color map
    if args.part == "real":
        to_be_plotted = chiqw_re
        label_chiqw = r"Re \chi(q,w)"
    elif args.part == "imag":
        to_be_plotted = chiqw_im
        label_chiqw = r"Im \chi(q,w)"
    else:
        raise ValueError(f"Invalid part: {args.part} (must be 'real' or 'imag')")

    if args.data_out is not None:
        with open(args.data_out, "w") as f:
            f.write("# x(q) omega Re[chi(q,omega)] Im[chi(q,omega)]\n")
            for iq in range(n_q):
                for iw in range(n_w):
                    f.write(f"{xarray[iq]} {ws[iw]} {chiqw_re[iq, iw]} {chiqw_im[iq, iw]}\n")
        print(f"Data written to {args.data_out}")

    im = ax.pcolormesh(X, Y, to_be_plotted.T, shading='auto', cmap='RdBu_r', norm=colors.CenteredNorm())
    fig.colorbar(im, ax=ax, label=label_chiqw)

    # Set x-axis ticks and labels
    ax.set_xticks(xticks)
    ax.set_xticklabels(xlabels)

    # Draw vertical lines at high symmetry points
    for x in xticks:
        ax.axvline(x=x, color='k', linestyle='-', alpha=0.3)

    ax.set_xlabel('q')
    ax.set_ylabel('ω')

    plt.tight_layout()
    for fmt in args.format:
        fmt = fmt.strip().lower()
        plt.savefig(f"{args.output_file_prefix}.{fmt}")
        print(f"Figure saved to {args.output_file_prefix}.{fmt}")


if __name__ == "__main__":
    main()
