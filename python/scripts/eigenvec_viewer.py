#!/usr/bin/env python3

import numpy as np
import os
import sys
import argparse
import math
import ast

from chiq.point_group import PointGroup


def print_matrix(mat):
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            z = mat[i,j]
            print("  %6.3lf" % z.real if abs(z.real) > 1e-4 else "   0    ", end="")
            print(" %+6.3lfj" % z.imag if abs(z.imag) > 1e-4 else "        ", end="")
        print("")


def print_vec_simply(v):
    def print_float(_x):
        # print "%.10f" %_x if abs(_x)>1e-9 else "0",
        if abs(_x)>1e-10:
            print(" {}".format(_x), end="")
        else:
            print(" 0", end="")
        # print(" {}".format(np.round(_x, 10)), end="")

    for z in v:
        print_float(z.real)
        print_float(z.imag)
        print(" ", end="")
    print("")


def print_weights(weights):
    for w in weights:
        if w == 0:
            # print("   .", end="")
            print("    .", end="")
        else:
            # print(" %3.0f" %(w*100.), end="")
            print("  %3.0f" %(w*100.), end="")
    print("")


# def find_max_rank(weights):
#     # for irp, w_irp in weights.items():
#     #     return len(w_irp)
#
#     max_rank = 0
#     for w_irp in weights.values():
#         for rank, w in enumerate(w_irp):
#             if w != 0:
#                 max_rank = max(max_rank, rank)
#     return max_rank

def find_max_rank(multipoles):
    max_rank = 0
    for m_irp in list(multipoles.values()):
        for m in m_irp:
            max_rank = max(max_rank, m.rank)
    return max_rank


def main():
    P = argparse.ArgumentParser()
    P.add_argument('filein', help="File name, e.g., chi_q_eigenvec.00.00.00.dat")
    # P.add_argument('-b', '--blocks', default=None, type=int, help="Number of blocks. If this is given, each block is printed as a matrix.")
    P.add_argument('-a', '--atoms', default=None, type=int, help="Number of atoms. If this is given, eigenvector for each atom is printed as a matrix.")
    # P.add_argument('-s', '--spin', action='store_true', help="Spin is conserved")
    group = P.add_mutually_exclusive_group()
    group.add_argument('--spin-orbit', action='store_true', help="Spin-orbit coupling is included in DMFT calculation. The spin is not conserved.")
    group.add_argument('--spin-charge', action='store_true', help="If this is activated, up-up and down-down components are converted to charge and S^z basis.")
    P.add_argument('-w', '--weight', choices=['total', 'rank'], default=None, help="print weight of each irrep in eigenvectors.")
    P.add_argument('-g', '--group', default=None, help="Group name (required if -w is activated)")
    P.add_argument('-l', '--l', default=None, type=str, help="Angular momentum, e.g., 2, 3, 5/2 (required if -w is activated)")
    P.add_argument('-t', '--time_reversal', action='store_true', help="Classify multipole operators by time-reversal symmetry (only when -w is activated)")
    P.add_argument('-f', '--fileout', action='store_true', help="Save weights into file (only when -w is activated)")
    P.add_argument('--basis_order', default=None, type=str, help="Order of basis, e.g., [2,1,0]")

    args = P.parse_args()
    print(args)

    # Load data file
    if not os.path.exists(args.filein):
        sys.exit(f"ERROR: File '{args.filein}' not found")

    print(f"\nLoad {args.filein!r}")
    row_vecs = np.loadtxt(args.filein, ndmin=2).view(complex)
    print(f" (n_eigen, dim) = {row_vecs.shape}")
    _, dim = row_vecs.shape

    # TODO: print information of multipoles

    print(f"\nDecompose dim={dim}")

    # -a option (atom)
    n_atoms = args.atoms
    if args.atoms is None:
        n_atoms = 1
    print(f" n_atoms = {n_atoms}")

    # --spin-orbit, --spin-charge
    flag_spin_orbit = args.spin_orbit
    flag_spin_charge = args.spin_charge
    if flag_spin_orbit:
        n_spin = 1
        str_ch_sp = [""]
    else:
        n_spin = 2
        if flag_spin_charge:
            str_ch_sp = ["ch", "sp+", "sp-", "spz"]  # charge/spin label
        else:
            str_ch_sp = ["uu", "ud", "du", "dd"]
    print(f" n_spin = {n_spin}")

    # n_orb
    n_orb = np.sqrt(dim / (n_atoms * n_spin**2))
    if n_orb.is_integer() is False:
        sys.exit(f"ERROR: can not determine n_orb by dim = n_atoms * (n_spin)^2 * (n_orb)^2.")
    n_orb = n_orb.round().astype(int)
    print(f" n_orb = {n_orb}")

    print(f" (atom, s1, s2, m1, m2) = ({n_atoms}, {n_spin}, {n_spin}, {n_orb}, {n_orb})")

    # -w option
    if args.weight is not None:
        if args.group is None or args.l is None:
            sys.exit("ERROR: -g (--group) and -l (--l) are required when -w (--weight) is given")

        # -t option
        time_reversal = args.time_reversal
        if time_reversal and args.weight != 'total':
            sys.exit("ERROR: -t (--time_reversal) option is valid only when --weight='total'")

        if args.basis_order is not None:
            args.basis_order = ast.literal_eval(args.basis_order)

        grp = PointGroup(args.l, args.group, basis_order=args.basis_order, verbose=False)
        basis = grp.basis
        print(f"\nPointGroup l={args.l}, group={args.group}")
        print(f" Bases: Jz = {basis}")

        # -f option
        def write(*args, **kwarge):  # dummy function
            pass
        if args.fileout:
            # fileout = os.path.splitext(args.filein)[0] + '_weight.dat'
            fileout = os.path.splitext(args.filein)[0] + '.weight'
            # fileout = args.filein.replace('eigenvec', 'weight')
            print(f" Weights are saved to {fileout!r}")
            f = open(fileout, 'w')
            write = f.write

            write(f"# PointGroup: l={args.l}, group={args.group!r}, basis={basis}\n")
            write(f"# n_atoms = {n_atoms}\n")
            # write(f"# n_spin = {n_spin}\n")
            write(f"# spin/charge = {str_ch_sp}\n")
            # write(f"# Bases: Jz = {basis}\n")
            write(f"# irreps = {grp.irreps1}\n")
            write(f"# {n_atoms*n_spin**2*len(grp.irreps1)} columns: Weights for (atom, spin/chage, irrep)\n")

    print("")
    if args.atoms is None:
        # print as row vector
        for i, vec in enumerate(row_vecs):
            print("#%d" %i)
            # print vec
            print_vec_simply(vec)
    else:
        n_atoms = args.atoms

        # print as matrix
        for i, vec in enumerate(row_vecs):
            print("=============")
            print("#%d" %i)

            # [a, s1, s2, m1, m2]
            vec_assmm = vec.reshape(n_atoms, n_spin, n_spin, n_orb, n_orb)

            # [a, ss, m1, m2]
            vec_asmm = vec_assmm.reshape(n_atoms, n_spin**2, n_orb, n_orb)

            if flag_spin_charge:
                # uu, ud, du, dd -> ch, sp+, sp-, spz
                for a in range(n_atoms):
                    ch = (vec_assmm[a, 0, 0, :, :] + vec_assmm[a, 1, 1, :, :]) / np.sqrt(2)
                    sp = (vec_assmm[a, 0, 0, :, :] - vec_assmm[a, 1, 1, :, :]) / np.sqrt(2)
                    vec_assmm[a, 0, 0, :, :] = ch
                    vec_assmm[a, 1, 1, :, :] = sp

            if args.weight is None:
                for a in range(n_atoms):
                    print(f"atom {a}")
                    for ss in range(n_spin**2):
                        print(f"{str_ch_sp[ss]}")
                        print_matrix(vec_asmm[a, ss, :, :])
                    print("")
            else:
                # grp = PointGroup(args.l, args.group, verbose=False)

                # print names of irreps
                print("             ", end="")
                if time_reversal:
                    for irp in grp.irreps1:
                        # print("%3s+" % irp, end="")
                        print(" %3s+" % irp, end="")
                    for irp in grp.irreps1:
                        # print("%3s-" % irp, end="")
                        print(" %3s-" % irp, end="")
                else:
                    for irp in grp.irreps1:
                        # print(" %3s" % irp, end="")
                        print("  %3s" % irp, end="")
                print("")
                # print("----------------------------------------")

                write(f"#{i}\n")

                # for each atom/spin
                for a in range(n_atoms):
                    for ss in range(n_spin**2):
                        mat = vec_asmm[a, ss, :, :]
                        weights_irps = grp.weights_mat(mat, time_reversal=time_reversal)
                        print(f" atom{a:2d} ", end="")
                        print(f"{str_ch_sp[ss]:3} |", end="")

                        print_weights(weights_irps)

                        for w in weights_irps:
                            write(f"{w} ")

                        if args.weight == 'total':
                            continue

                        # print rank-resolved weight
                        print("             ", end="")
                        for irp in grp.irreps1:
                            print("----", end="")
                        print()

                        weights = grp.weights_mat_multipoles(mat)
                        max_rank = find_max_rank(grp.multipoles)
                        for rank in range(max_rank+1):
                            print("     rank%2d |" % rank, end="")
                            weights_irps = [weights[irp][rank] for irp in grp.irreps1]
                            print_weights(weights_irps)
                        print("")
                    write("\n")

    if 'f' in vars():  # f is defined
        f.close()


if __name__ == '__main__':
    main()
