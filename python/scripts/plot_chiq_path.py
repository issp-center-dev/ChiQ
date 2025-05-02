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


def marker(i):
    markers = ["-o", "-v", "-^", "-*", "-<", "->", "-D", "-d", "-h", "-p"]
    return markers[(i // 10) % len(markers)]


_cmap = plt.get_cmap("tab10")  # default color map


def cmap(i):
    return _cmap(i % _cmap.N)


def main():
    P = argparse.ArgumentParser()
    P.add_argument("file_qpath")
    P.add_argument("file_eigen")
    P.add_argument(
        "-d", action="store_true", help="replace \chi with \Delta\chi in y-label"
    )
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
    # P.add_argument('-i', action='store_true', help="plot inverse susceptibility")
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
    E = ChiQEigenPath(args.file_qpath, wb=args.wb)
    xarray = E.get_x()
    yarray = E.get_y_on_path(args.file_eigen)

    assert isinstance(xarray, np.ndarray)
    assert isinstance(yarray, np.ndarray)
    n_q, n_col = yarray.shape
    print("# of q points", n_q)
    print("# of columns", n_col)
    assert xarray.shape == (n_q,)
    assert yarray.shape == (n_q, n_col)

    # --------------------------------------------------------------------------
    # mode dependent variables
    dict_flag_inv = {
        "chi": True,
        "chi0": False,
        "rpa": True,
        "scl": True,
        "rrpa": True,
        "Iq": False,
    }
    # dict_prefix = {'chi': 'chi_q', 'chi0': 'chi0_q', 'rpa': 'chi_q_rpa', 'scl': 'chi_q_scl', 'rrpa': 'chi_q_rrpa', 'Iq': 'I_q'}

    flag_inv = dict_flag_inv[args.mode]
    # prefix = dict_prefix[args.mode]

    filein_base = os.path.splitext(args.file_eigen)[0]

    # --------------------------------------------------------------------------
    # save data
    file_data = args.data_out
    if file_data is not None:
        xticks, xlabels = E.get_xticks(latex=False)

        # information of symmetric point
        header = ""
        for _x, _label in zip(xticks, xlabels):
            header += f'"{_label}" {str(_x)}, '

        # xarray [N_k]
        # yarray [N_k, N_mode]
        np.savetxt(file_data, np.hstack((xarray[:, None], yarray)), header=header)
        print("'%s'" % file_data)

    # --------------------------------------------------------------------------
    # labels
    mode_labels = {}
    if args.label is None:
        mode_labels = {i: str(i) for i in range(n_col)}
    else:
        with open(args.label, "r") as f:
            txt = f.read()
            # print(r"%s" %txt)
            mode_labels = ast.literal_eval(txt)
            assert isinstance(mode_labels, dict)
        print("mode_labels =")
        for key, label in list(mode_labels.items()):
            print(f"  {key} : {label}")

    # --------------------------------------------------------------------------
    # configuration of subfigures
    cols_list = []
    if args.subfigures is not None:
        with open(args.subfigures, "r") as f:
            for line in f:
                if line[0] == "#":  # skip comment
                    continue
                array = line.split()
                if len(array) == 0:  # skip empty line
                    continue
                cols_list.append(np.array(array, dtype=int))
        print("cols_list =")
        for j, cols in enumerate(cols_list):
            print(f"  {j} : {cols}")

    # --------------------------------------------------------------------------
    # plot

    xticks, xlabels = E.get_xticks(latex=True)

    str_chi = r"\chi(q)"
    bgcolor = "lightyellow"
    if args.d:
        str_chi = r"\Delta" + str_chi
    if args.mode == "Iq":
        str_chi = r"I(q)"
        bgcolor = "oldlace"

    def plot_common(_ax, _fig, _filefig, _bgcolor="none"):
        # x ticks
        _ax.set_xticks(xticks)
        _ax.set_xticklabels(xlabels)
        _ax.grid(axis="x")
        _ax.set_axisbelow(True)
        _ax.margins(x=0)
        _ax.axhline(y=0, color="k", zorder=-1)  # show yzero
        _ax.legend(
            prop={"size": args.label_fontsize},
            bbox_to_anchor=(1.04, 1),
            loc="upper left",
        )
        _ax.set_facecolor(_bgcolor)

        # save
        _fig.tight_layout()
        _fig.savefig(_filefig)
        print(f"{_filefig!r}")

    # plot all
    fig, ax_all = plt.subplots()
    for i in range(n_col):
        # skip column i if label is not given
        if i in mode_labels:
            ax_all.plot(
                xarray,
                yarray[:, i],
                marker(i),
                color=cmap(i),
                label=mode_labels[i],
                zorder=5,
            )
    ax_all.set_ylabel(rf"${str_chi}$")
    ax_all.set_ylim(args.ymin, args.ymax)
    for fmt in args.format:
        plot_common(ax_all, fig, f"{filein_base}_path.{fmt}", _bgcolor=bgcolor)

    # plot 1/chi
    if flag_inv:
        fig, ax = plt.subplots()
        for i in range(n_col):
            if i in mode_labels:
                ax.plot(
                    xarray,
                    1.0 / yarray[:, i],
                    marker(i),
                    color=cmap(i),
                    label=mode_labels[i],
                )
        ax.set_ylabel(rf"$1 / {str_chi}$")
        ax.set_ylim(-0.4, 1)
        for fmt in args.format:
            plot_common(ax, fig, f"{filein_base}_path_inv.{fmt}", _bgcolor="lightgray")

    # plot parts (optional)
    for j, cols in enumerate(cols_list):
        fig, ax = plt.subplots()
        for i in cols:
            # skip column i if label is not given
            if i in mode_labels:
                ax.plot(
                    xarray,
                    yarray[:, i],
                    marker(i),
                    color=cmap(i),
                    label=mode_labels[i],
                    zorder=5,
                )
        ax.set_ylabel(rf"${str_chi}$")
        ax.set_ylim(args.ymin, args.ymax)
        if args.sharey:
            ax.sharey(ax_all)
        for fmt in args.format:
            plot_common(ax, fig, f"{filein_base}_path_{j}.{fmt}", _bgcolor=bgcolor)


if __name__ == "__main__":
    main()
