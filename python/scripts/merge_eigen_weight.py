#!/usr/bin/env python3

import os
import sys
import argparse
import numpy as np


def _read_weights(_file):
    weights_list = []
    with open(_file, 'r') as f:
        for line in f:
            if line[0] != '#':  # skip comment line
                weights_list.append(line)
    return weights_list


def main():
    P = argparse.ArgumentParser()
    P.add_argument('filein', help="File name, e.g., 'chi_q_eigen.dat'")
    # P.add_argument('-q', '--q_path', type=str, help="q-path file, e.g., 'q_path.dat'. If this is given, symbols like '00.00.00' are replaced with x-coordinate for q-path plot.")

    args = P.parse_args()
    print(args)

    # Load data file
    filein = args.filein
    if not os.path.exists(filein):
        sys.exit(f"ERROR: File '{filein}' not found")
    print(f"\nfilein = {filein!r}")

    # output filename
    filein_base = os.path.splitext(filein)[0]
    fileout = f"{filein_base}.weight"
    print(f"fileout = {fileout!r}")

    # basename of input weight file
    fileweight_base = os.path.splitext(filein)[0] + 'vec'
    # fileweight_base = os.path.splitext(filein.replace('eigen', 'weight'))[0]
    # print(fileweight_base)

    # get n_eig
    data_filein = np.loadtxt(filein, dtype=str)
    n_q, n_eig = data_filein.shape
    n_eig -= 2  # minus 2 columns (w, q)
    del data_filein
    print(f"# of eigenvalues : {n_eig}")
    print(f"# of q-points : {n_q}")

    print("")
    with open(filein, "r") as fin, open(fileout, "w") as fout:
        for line in fin:
            if line[0] == '#':  # skip comment line
                continue
            array = line.split()
            # print(array)
            w = array[0]
            q = array[1]
            eigenvalues = array[2:]
            print("q =", q)
            assert len(eigenvalues) == n_eig

            # fileweight = f"{fileweight_base}.{q}_weight.dat"
            fileweight = f"{fileweight_base}.{q}.weight"
            print(f"  {fileweight}", end="")

            if os.path.exists(fileweight):
                print(f" -> found")

                weights = _read_weights(fileweight)
                assert len(weights) == n_eig

                for i, eigv in enumerate(eigenvalues):
                    fout.write(f"{w} {q} {eigv} {weights[i]}")  # weights contain line break
            else:
                print(" -> not found")
                continue
                for _, eigv in enumerate(eigenvalues):
                    fout.write(f"{w} {q} {eigv}\n")  # weights contain line break
            fout.write("\n")


if __name__ == '__main__':
    main()
