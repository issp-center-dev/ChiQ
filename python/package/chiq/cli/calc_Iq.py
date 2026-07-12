#!/usr/bin/env python3

import numpy as np
from scipy import linalg
import sys
import argparse
from itertools import product

from chiq.h5bse import h5BSE
from chiq.matrix_dict import MatrixDict


def _is_zero(mat):
    return np.all(np.absolute(mat) < 1e-8)


def _convert_to_matrix(block_matrix, n_block, n_inner):
    assert isinstance(block_matrix, dict)
    dim = n_inner * n_block
    mat = np.zeros((n_block, n_inner, n_block, n_inner), dtype=complex)
    for block, bmat in block_matrix.items():
        i, j = block
        mat[i, :, j, :] = bmat
    return mat.reshape((dim, dim))


def _convert_to_dict(mat, n_block, n_inner):
    mat4 = mat.reshape((n_block, n_inner, n_block, n_inner))
    block_matrix = {}
    for i, j in product(range(n_block), repeat=2):
        # print(i, j)
        bmat = mat4[i, :, j, :]
        assert bmat.shape == (n_inner, n_inner)

        if not _is_zero(bmat):
            block_matrix[(i, j)] = bmat

    return block_matrix


def _eigh_descending(mat):
    eigvals, eigvecs = linalg.eigh(mat)
    return eigvals[::-1], eigvecs[:, ::-1]


def _is_hermitian(mat):
    antihermite = (mat - mat.conjugate().T) / 2
    return _is_zero(antihermite)


def _hermitianize(mat):
    return (mat + mat.conjugate().T) / 2


def _load_chi_loc(BS, w):
    chi_loc = BS.get(key=('chi_loc', w))
    if chi_loc is False:
        print("  'chi_loc' not found. Compute 'chi_loc' from 'X_loc'.")
        # Compute chi_loc from X_loc
        beta = BS.get(key=('beta',))
        x_loc = BS.get(key=('X_loc', w))
        if x_loc is False:
            sys.exit("Either chi_loc or X_loc is required.")
        chi_loc = {key: np.sum(array, axis=(2, 3))/beta for key, array in list(x_loc.items())}
    else:
        print("  'chi_loc' found")
    return chi_loc


def calc_Iq(h5_file : str, wb : np.ndarray, retain : int, remove : int, verbose=True):

    print(f"\nLoad data from {h5_file!r}")
    BS = h5BSE(h5_file)
    BS.open('r')  # Keep HDF5 file open to improve performance. Close manually.

    inner_name = BS.get(key=('inner_name'))
    block_name = BS.get(key=('block_name'))
    print(f"  inner_name: {inner_name}")
    print(f"  block_name: {block_name}")
    n_inner = len(inner_name)
    n_block = len(block_name)
    dim = n_inner * n_block

    # -------------------------------------------------------------
    # chi_loc

    print(f"\nLoad 'chi_loc'")
    chi_loc = _load_chi_loc(BS, wb)
    assert isinstance(chi_loc, dict)

    # Convert block matrix (dict) to ndarray
    mat_chi_loc = _convert_to_matrix(chi_loc, n_block, n_inner)
    assert mat_chi_loc.shape == (dim, dim)

    # hermitianize chi_loc
    if not _is_hermitian(mat_chi_loc):
        print("Warning: chi_loc is not hermitian.", file=sys.stderr)
    mat_chi_loc = _hermitianize(mat_chi_loc)

    # Get eigenvalues in descending order
    print("\nDiagonalize chi_loc")
    eigs, u = _eigh_descending(mat_chi_loc)
    assert eigs.shape == (dim,)
    assert u.shape == (dim, dim)

    # Number of eigenmodes retained
    n_retain = retain
    if n_retain is None:
        n_retain = dim - remove
    print(f"  Number of eigenmodes retained = {n_retain}")

    for i, eig in enumerate(eigs):
        if i == n_retain:
            print(f"  ------------------------- Remove eigenmodes below")
        print(f"  {i}: {eig}")

    # Dimensionality reduction
    eigs = eigs[:n_retain]
    u = u[:, :n_retain]
    u_dag = u.conj().T
    assert eigs.shape == (n_retain,)

    # -------------------------------------------------------------
    # chi_q

    # Get list of 'chi_q' data
    keylist = BS.get_keylist_data(input_output='output')
    keylist_chiq = [key for key in keylist if key[0:2]==('chi_q', wb)]
    BS.close()

    print(f"\n{len(keylist_chiq)} data of 'chi_q' found")
    print(f"Compute 'I_q' (existing data are overwritten)")

    print(f"\nStart q-loop")
    BS.open('a')  # Keep HDF5 file open to improve performance. Close manually.
    for key in keylist_chiq:
        if verbose:
            print(f"Load data: key={key}")
        chi_q = BS.get(key=key)

        # Convert block matrix (dict) to ndarray
        mat_chi_q = _convert_to_matrix(chi_q, n_block, n_inner)

        # Compute I_q
        I_q = u @ (np.diag(eigs**(-1)) - linalg.inv(u_dag @ mat_chi_q @ u)) @ u_dag
        assert I_q.shape == (dim, dim)

        # Convert ndarray to block matrix (dict)
        I_q_dict = _convert_to_dict(I_q, n_block, n_inner)
        assert isinstance(I_q_dict, dict)

        _, _w, q = key
        assert _w == wb
        if verbose:
            print(f"Save data: key=('I_q', {wb}, {q})")
        BS.save(key=('I_q', wb, q), data=I_q_dict)
    BS.close()
    print("End q-loop")


def main():
    P = argparse.ArgumentParser()
    P.add_argument('-f', '--file', default='dmft_bse.out.h5', help="h5 filename. Default is 'dmft_bse.out.h5'")
    group = P.add_mutually_exclusive_group()
    group.add_argument('--remove', type=int, default=1, help="Number of eigenmodes removed. Default is 1, which intends to remove the charge susceptibility.")
    group.add_argument('--retain', type=int, default=None, help="Number of eigenmodes retained.")
    P.add_argument('-w', type=int, default=0, help="Bosonic Matsubara frequency. Default is 0 (static susceptibility).")
    P.add_argument('--verbose', action='store_true', help="verbose output")
    args = P.parse_args()
    print(args)

    calc_Iq(h5_file=args.file, wb=args.w, retain=args.retain, remove=args.remove, verbose=args.verbose)

if __name__ == "__main__":
    main()
