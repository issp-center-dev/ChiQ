#!/usr/bin/env python3

import numpy as np
import os
import sys
import argparse
import ast
from collections import namedtuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def marker(i):
    markers = ["o", "v", "^", "*", "<", ">", "D", "d", "h", "p"]
    return markers[(i // 10) % len(markers)]


_cmap = plt.get_cmap("tab10")  # default color map


def cmap(i):
    return _cmap(i % _cmap.N)


class IrEigen4Plot(object):
    """
    Generate list of (x,y) coodinates to be plotted

    How to use:
        E = IrEigen4Plot(file_eigen)
        xarray = E.get_x()
        yarray = E.get_y()
    """

    def __init__(self, _file_eigen):
        self.__set_suscep(_file_eigen)

    def __set_suscep(self, _file_eigen):
        """Read I(q) data"""
        Iq = namedtuple("Iq", ["w", "q", "eigen"])

        self._data = []  # suscep[(w,q)] = [suscep1, suscep2, ...]
        for array in np.loadtxt(_file_eigen, dtype=str, comments="#", ndmin=2):
            w = int(array[0])
            qvec = np.array(array[1].split("."), dtype=int)
            eigen = np.array(array[2:], dtype=float)
            self._data.append(Iq(w, qvec, eigen))

    def get_x(self, w, avec=None):
        """
        Get array of radial coordinate r for plot

        w: int, bosonic frequency
        avec: np.ndarray(3,3), lattice vectors (row vectors)

        Return
            x: np.ndarray(n_r)
        """
        if avec is None:
            avec = np.identity(3)
        assert isinstance(avec, np.ndarray)
        assert avec.shape == (3, 3)

        x = []
        for iq in self._data:
            if iq.w == w:
                # x.append(np.linalg.norm(iq.q))
                x.append(np.linalg.norm(iq.q @ avec))  # q_i * A_{ij}
        return np.array(x)

    def get_y(self, w):
        """
        Get array of y values

        Return
            y: np.ndarray(n_r, n_mode)
        """
        y = []
        for iq in self._data:
            if iq.w == w:
                y.append(iq.eigen)
        return np.array(y)


# =============================================================================


def main():

    P = argparse.ArgumentParser()
    P.add_argument("file_eigen")
    # P.add_argument('--mode', default='I_r', choices=['I_r',], help="default: I_r")
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
    P.add_argument(
        "--avec",
        default=None,
        help="set filename to get the lattice vectors. The file should contain a real 3x3 matrix with row vectors a1, a2, a3. Default is a unit matrix.",
    )
    P.add_argument(
        "--subfigures",
        "--subfigure",
        default=None,
        help="set filename to define optional figures that show a part of data. One line in the file defines one figure. Each line should contain space-separated integers, which specify columns plotted in the figure.",
    )
    P.add_argument("--xmin", default=0, type=float, help="Left bound of x axis")
    P.add_argument("--xmax", default=5, type=float, help="Right bound of x axis")
    P.add_argument("--ymin", default=None, type=float, help="Upper bound of y axis")
    P.add_argument("--ymax", default=None, type=float, help="Lower bound of y axis")
    P.add_argument("-w", default=0, type=int, help="bosonic frequency")
    P.add_argument(
        "--sharey",
        action="store_true",
        help="share y-axis (all figures have the same y range)",
    )
    args = P.parse_args()

    print("\nRunning", os.path.basename(__file__))
    print(args)

    # --------------------------------------------------------------------------
    # lattice vectors
    if args.avec is None:
        avec = np.identity(3)
    else:
        try:
            print(f"Loading avec from '{args.avec}'")
            avec = np.loadtxt(args.avec)
            assert avec.shape == (3, 3), f"avec.shape={avec.shape} must be (3, 3)"
        except Exception as e:
            print("Error:", e, file=sys.stderr)
            exit(1)
    print(f"avec = {avec.tolist()}")  # print in one line

    # --------------------------------------------------------------------------
    # get data to plot
    E = IrEigen4Plot(args.file_eigen)
    xarray = E.get_x(args.w, avec)
    yarray = E.get_y(args.w)

    assert isinstance(xarray, np.ndarray)
    assert isinstance(yarray, np.ndarray)
    n_r, n_col = yarray.shape
    print("# of r points", n_r)
    print("# of columns", n_col)
    assert xarray.shape == (n_r,)
    assert yarray.shape == (n_r, n_col)

    # --------------------------------------------------------------------------
    # mode dependent variables
    # dict_prefix = {'I_r': 'I_r'}

    # prefix = dict_prefix[args.mode]
    filein_base = os.path.splitext(args.file_eigen)[0]

    # --------------------------------------------------------------------------
    # save data
    file_data = args.data_out
    if file_data is not None:
        # xarray [N_k]
        # yarray [N_k, N_mode]
        np.savetxt(file_data, np.hstack((xarray[:, None], yarray)))
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
            print(" ", key, ":", label)

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

    def plot_common(_ax, _fig, _filefig, _bgcolor="none"):
        # x ticks
        # _ax.set_xticks(xticks)
        # _ax.set_xticklabels(xlabels)
        _ax.grid(axis="x")
        _ax.set_xlabel(r"$|r|$")
        _ax.set_ylabel(r"$I(r)$")
        _ax.set_xlim(args.xmin, args.xmax)
        _ax.set_ylim(args.ymin, args.ymax)
        _ax.set_axisbelow(True)
        # _ax.margins(x=0)
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

    bgcolor = "lightcyan"

    # plot I(r)
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
    for fmt in args.format:
        plot_common(ax_all, fig, f"{filein_base}_distance.{fmt}", _bgcolor=bgcolor)

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
        if args.sharey:
            ax.sharey(ax_all)
        for fmt in args.format:
            plot_common(ax, fig, f"{filein_base}_distance_{j}.{fmt}", _bgcolor=bgcolor)


if __name__ == "__main__":
    main()
