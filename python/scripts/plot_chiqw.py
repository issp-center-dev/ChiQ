#!/usr/bin/env python3

import numpy as np
import sys
import os
import argparse
import ast

from chiq.chiqw_eigen_path import ChiQWEigenPath

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def marker(i):
    markers = ["-o", "-v", "-^", "-*", "-<", "->", "-D", "-d", "-h", "-p"]
    return markers[(i // 10) % len(markers)]


_cmap = plt.get_cmap("tab10")  # default color map


def cmap(i):
    return _cmap(i % _cmap.N)


def main():
    P = argparse.ArgumentParser()
    P.add_argument("file_qpath")
    P.add_argument("dir_chiqw")
    P.add_argument("-f", "--file_prefix", default="chiq_w")
    P.add_argument("-i", "--ielem", type=int, default=0)
    P.add_argument("--wmin", type=float, default=None)
    P.add_argument("--wmax", type=float, default=None)
    P.add_argument(
        "--data_out", default=None, help="set filename to get numerical data for plot"
    )
    P.add_argument(
        "--format",
        default="pdf",
        type=lambda s: s.split(","),
        help="output format. multiple formats can be specified by comma-separated string, e.g., 'pdf,png'.",
    )
    # P.add_argument('-f', '--file', const="suscep", nargs="?", help="Base filename of pdf")
    P.add_argument(
        "--label",
        "-l",
        default=None,
        help="set filename to get labels for chi(q) plot. The file should contain {0: 'label0', 1: 'label1'} in dict format. Skipped entries will not shown in plot.",
    )
    P.add_argument(
        "--label-fontsize", default=8, type=int, help="set fontsize of labels"
    )
    P.add_argument(
        "--subfigures",
        "--subfigure",
        default=None,
        help="set filename to define optional figures that show a part of data. One line in the file defines one figure. Each line should contain space-separated integers, which specify columns plotted in the figure.",
    )
    P.add_argument("--ymin", default=None, type=float, help="Lower bound of y axis")
    P.add_argument("--ymax", default=None, type=float, help="Upper bound of y axis")
    P.add_argument(
        "--sharey",
        action="store_true",
        help="share y-axis (all figures have the same y range)",
    )
    P.add_argument("-w", "--wb", default=0, type=int, help="Index of Matsubara frequency to be plotted")
    args = P.parse_args()

    print("\nRunning", os.path.basename(__file__))
    print(args)

    # --------------------------------------------------------------------------
    # get data to plot
    E = ChiQWEigenPath(args.file_qpath)
    xarray = E.get_x()
    ws, chiqw_re, chiqw_im = E.get_chiqw_on_path(dir_chiqw=args.dir_chiqw,
                                             file_prefix=args.file_prefix,
                                             ielem=args.ielem, wmin=args.wmin,
                                             wmax=args.wmax)

    n_q, n_w = chiqw_re.shape
    print("# of q points", n_q)
    print("# of w points", n_w)
    assert xarray.shape == (n_q,)
    assert chiqw_re.shape == (n_q, n_w)
    assert chiqw_im.shape == (n_q, n_w)

    # --------------------------------------------------------------------------
    # plot

    xticks, xlabels = E.get_xticks(latex=True)

    str_chiqw = r"\chi(q,w)"

    fig, ax = plt.subplots(figsize=(8, 6))

    # Create meshgrid for pcolormesh
    X, Y = np.meshgrid(xarray, ws)
    
    # Plot color map
    im = ax.pcolormesh(X, Y, chiqw_re.T, shading='auto', cmap='RdBu_r')
    fig.colorbar(im, ax=ax, label=f'Re {str_chiqw}')

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
        plt.savefig(f"chiqw.{fmt}")


if __name__ == "__main__":
    main()
