#!/usr/bin/env python3

import numpy as np
import sys
import os
import argparse
import ast

from chiq.chiq_eigen_path import ChiQEigenPath

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


_cmap = plt.get_cmap("tab10")  # default color map


def cmap(i):
    return _cmap(i % _cmap.N)


def main():
    P = argparse.ArgumentParser()
    P.add_argument("file_qpath")
    P.add_argument("file_eigen_weight")
    P.add_argument(
        "--mode",
        default="chi",
        choices=["chi", "chi0", "rpa", "scl", "rrpa", "Iq"],
        help="default: chi",
    )
    P.add_argument(
        "--data_out", default=None, help="set filename to get numerical data for plot"
    )
    P.add_argument(
        "--format",
        default="pdf",
        type=lambda s: s.split(","),
        help="output format. multiple formats can be specified by comma-separated string, e.g., 'pdf,png'.",
    )
    P.add_argument(
        "--label",
        "-l",
        default=None,
        help="set filename to get labels for chi(q) plot. The file should contain {0: 'label0', 1: 'label1'} in dict format. Skipped entries will not shown in plot.",
    )
    P.add_argument(
        "--label-fontsize", default=8, type=int, help="set fontsize of labels"
    )
    P.add_argument("--ymin", default=None, type=float, help="Lower bound of y axis")
    P.add_argument("--ymax", default=None, type=float, help="Upper bound of y axis")
    args = P.parse_args()

    print("\nRunning", os.path.basename(__file__))
    print(args)

    # --------------------------------------------------------------------------
    # get data to plot
    E = ChiQEigenPath(args.file_qpath)
    xarray = E.get_x()
    assert isinstance(xarray, np.ndarray)
    print(f"# of q points: {xarray.size}")

    points = E.get_points(args.file_eigen_weight)
    assert isinstance(points, np.ndarray)
    # print(f"points.shape = {points.shape}")

    # data for plot
    x = points[:, 0]
    y = points[:, 1]
    weights = points[:, 2:]
    n_weights = weights.shape[1]
    print(f"# of points to plot: {points.shape[0]}")
    print(f"# of weights: {n_weights}")

    # --------------------------------------------------------------------------
    # output filename
    filein_base = os.path.splitext(args.file_eigen_weight)[0]

    # --------------------------------------------------------------------------
    # save data
    file_data = args.data_out
    if file_data is not None:
        xticks, xlabels = E.get_xticks(latex=False)

        # information of symmetric point
        header = ""
        for _x, _label in zip(xticks, xlabels):
            header += f'"{_label}" {str(_x)}, '

        np.savetxt(file_data, points, header=header)
        print("'%s'" % file_data)

    # --------------------------------------------------------------------------
    # labels
    mode_labels = {}
    if args.label is None:
        mode_labels = {i: str(i) for i in range(n_weights)}
    else:
        with open(args.label, "r") as f:
            txt = f.read()
            # print(r"%s" %txt)
            mode_labels = ast.literal_eval(txt)
            assert isinstance(mode_labels, dict)
        print("mode_labels =")
        for key, label in list(mode_labels.items()):
            print(" ", key, ":", label)

    # --------------------------------------------------------------------------
    # plot

    xticks, xlabels = E.get_xticks(latex=True)

    ylabel = r"$\chi(q)$"
    if args.mode == "Iq":
        ylabel = r"$I(q)$"

    def plot_common(_ax, _fig, _filefig):
        # x ticks
        _ax.set_xticks(xticks)
        _ax.set_xticklabels(xlabels)
        _ax.grid(axis="x")
        _ax.set_axisbelow(True)
        _ax.margins(x=0)
        _ax.set_ylabel(ylabel)
        _ax.axhline(y=0, color="k", zorder=-1)  # show yzero
        _ax.legend(
            prop={"size": args.label_fontsize},
            bbox_to_anchor=(1.04, 1),
            loc="upper left",
        )  # updated fontsize

        # save
        _fig.tight_layout()
        _fig.savefig(_filefig)
        print(f"{_filefig!r}")

    opt = dict(
        # clip_on = False,
        zorder=5,
        linewidths=(1,),
        # edgecolors = ('black',),
        alpha=0.8,
    )

    # Plot all
    fig, ax_all = plt.subplots()
    for j in range(n_weights):
        ax_all.scatter(
            x=x,
            y=y,
            s=50 * weights[:, j] ** 2,
            label=mode_labels[j],
            color=cmap(j),
            **opt,
        )
    ax_all.set_ylim(args.ymin, args.ymax)
    for fmt in args.format:
        plot_common(ax_all, fig, f"{filein_base}_weight.{fmt}")

    # Plot each
    for j in range(n_weights):
        fig, ax = plt.subplots()
        ax.sharey(ax_all)
        ax.scatter(
            x=x,
            y=y,
            s=50 * weights[:, j] ** 2,
            label=mode_labels[j],
            color=cmap(j),
            **opt,
        )
        for fmt in args.format:
            plot_common(ax, fig, f"{filein_base}_weight_{j}.{fmt}")


if __name__ == "__main__":
    main()
