#!/usr/bin/env python3

import numpy as np
import argparse
from collections import Counter
import ast
import sys
import os

from chiq.point_group import PointGroup


def print_cf_states(_cf_states):
    print("\ncrystal field states")
    print(" name  dim  shape_of_vecs")
    for cf in _cf_states:
        # print("   {:<3}   {:>3}  {:>10}".format(cf.name, cf.dim, cf.vecs.shape))
        print("   {}   {}  {}".format(cf.name, cf.dim, cf.vecs.shape))


def output_cf_states(_filename, _cf_states):
    print("\noutput cf_states...")
    with open(_filename, 'w') as f:
        for cf in _cf_states:
            print("# irp={}, dim={}".format(cf.name, cf.dim), file=f)
            for col in range(cf.dim):
                print(np.round(cf.vecs[:, col], 10), file=f)
    print("'{}' created".format(_filename))


def print_multipole_stat(_multipoles):
    print("\nmultipole operators")
    print(" irreps  num_of_ops")
    for irp, m_irp in list(_multipoles.items()):
        # print("   {:<3}   {:>3}".format(irrep, len(multipoles[irrep])))
        count_rank = Counter([m.rank for m in m_irp])
        # print(count_rank)
        if count_rank:
            print("%3s" % irp, end="")
            for rank, count in list(count_rank.items()):
                print("   (rank {}) * {}".format(rank, count), end="")
            print("")


def output_multipole_bases_for_bse(_filename, _multipoles, _add_spin):
    def expand_with_spin(array, which):
        zeros = np.zeros(array.shape, dtype=complex)
        sqrt2 = np.sqrt(2.)
        if which == '':
            return array
        if which == 'charge':
            return np.hstack((array / sqrt2, zeros, zeros, array / sqrt2))
        elif which == 'spin_z':
            return np.hstack((array / sqrt2, zeros, zeros, -array / sqrt2))
        elif which == 'spin_+':
            return np.hstack((zeros, array, zeros, zeros))
        elif which == 'spin_-':
            return np.hstack((zeros, zeros, array, zeros))
        else:
            raise Exception

    spin_components = ['charge', 'spin_z', 'spin_+', 'spin_-'] if _add_spin else ['']
    with open(_filename, 'w') as f:
        for ch_sp in spin_components:
            print("### {}".format(ch_sp), file=f)
            for m_irp in list(_multipoles.values()):
                for m in m_irp:
                    print("# {}".format(m), file=f)
                    m_1d = np.array(m.mat).flatten()
                    for z in expand_with_spin(m_1d, ch_sp):
                        print("{} {} ".format(z.real, z.imag), file=f, end="")
                    print("", file=f)
    print("'{}' created".format(_filename))


def output_multipole_labels_for_bse(_filename, _multipoles, _add_spin):
    # Output in dictionary format, e.g.
    # {
    #   0: 'charge A1-rank0',
    #   1: 'SzSz A1--rank0',
    # }

    spin_components = ['charge ', 'SzSz ', 'S+S- ', 'S-S+ '] if _add_spin else ['']
    with open(_filename, 'w') as f:
        n = 0
        print("{", file=f)
        for ch_sp in spin_components:
            for m_irp in _multipoles.values():
                for m in m_irp:
                    # print(f"{n} '{ch_sp}{m.irp}-rank{m.rank}'", file=f)
                    print(f"  {n}: '{ch_sp}{m.irp}-rank{m.rank}',", file=f)
                    n += 1
        print("}", file=f)
    print("'{}' created".format(_filename))


# def main(angular_momentum, group, bases, fileout=False, filename_cf="cf_states.dat", filename_mp="multipole_bases.dat", basis_order=None, verbose=False):

def main():
    P = argparse.ArgumentParser(
        description="Print info of irreducible representations of crystal-field states and multipoles for a given set of (J, G), where J is an angular momentum and G is a point group. If '-f' option is activated, eigenbases are stored in files.",
        epilog="Available group_name : {}".format(PointGroup.group_list())
    )

    # P.add_argument('-l', '--list', action='store_true', help="Show the list of available group names")
    # args, _ = P.parse_known_args()
    # if args.list:
    #     print("Group list :", PointGroup.group_list())
    #     sys.exit(0)

    P.add_argument('angular_momentum', type=str, help="Angular momentum. Ex. 2, 3, 5/2.")
    P.add_argument('group_name', help="Name of point group. See below for the list of groups.")
    # P.add_argument('--filename', default="multipole_bases.dat", help="output filename")
    P.add_argument('-f', action='store_true', help="Data files are created")
    P.add_argument('-v', '--verbose', action='store_true')
    # group = P.add_mutually_exclusive_group()
    P.add_argument('--bases', default=None, help="List of irrep of bases. Ex. --bases=\"['E',]\"")
    P.add_argument('--basis_order', default=None, type=str, help="Order of basis. Ex. [2,1,0]")

    args = P.parse_args()
    print(args)

    fileout = args.f
    bases = args.bases
    filename_cf = "cf_states.dat"
    filename_mp = "multipole_bases.dat"
    filename_mp_original_basis = "multipole_bases_in_original_space.dat"
    filename_labels = "multipole_labels.dat"

    if args.basis_order is not None:
        args.basis_order = ast.literal_eval(args.basis_order)

    print("\n============= J = {} ================================".format(args.angular_momentum))
    grp = PointGroup(
        angular_momentum = args.angular_momentum,
        group = args.group_name,
        basis_order = args.basis_order,
        verbose = args.verbose,
    )

    print("irreps =", grp.irreps)

    # CEF states
    cf_states = grp.cf_states
    print_cf_states(cf_states)

    if fileout:
        output_cf_states(filename_cf, cf_states)

    if bases is not None:
        print("\nreduce space")
        list_bases = ast.literal_eval(bases)
        print(" {}".format(list_bases))
        cf_dict = {cf.name: cf for cf in cf_states}
        assert len(cf_states) == len(cf_dict)
        cfs = [cf_dict[irp] for irp in list_bases]

        # print(cfs)
        grp.reduce_states(cfs)

    # Multipoles
    print_multipole_stat(grp.multipoles)

    if fileout:
        add_spin = not grp.flag_half_integer  # True if angular_momentum is integer

        print("\noutput multipole labels for BSE...")
        output_multipole_labels_for_bse(filename_labels, grp.multipoles, add_spin)

        print("\noutput multipole bases for BSE...")
        output_multipole_bases_for_bse(filename_mp, grp.multipoles, add_spin)

        if bases is not None:
            # Output multipole operators represented in the original space
            output_multipole_bases_for_bse(filename_mp_original_basis, grp.multipoles_original_space, add_spin)


if __name__ == '__main__':
    main()
