#!/usr/bin/env python3

import numpy as np
from scipy import linalg
import sys
import os
import argparse
import configparser
from itertools import product
from collections import OrderedDict

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


def _convert_to_matrix_iw(block_matrix, n_block, n_inner, n_iw):
    assert isinstance(block_matrix, dict)
    dim = n_inner * n_block
    mat = np.zeros((n_iw, n_block, n_inner, n_block, n_inner), dtype=complex)
    for block, bmat in block_matrix.items():
        assert bmat.shape == (n_inner, n_inner, n_iw)
        i, j = block
        # [o1, o2, iw] -> [iw, o1, o2]
        mat[:, i, :, j, :] = bmat.transpose((2, 0, 1))
    return mat.reshape((n_iw, dim, dim))


def _convert_to_dict(mat, n_block, n_inner):
    mat4 = mat.reshape((n_block, n_inner, n_block, n_inner))
    block_matrix = {}
    for i, j in product(range(n_block), repeat=2):
        bmat = mat4[i, :, j, :]
        assert bmat.shape == (n_inner, n_inner)

        if not _is_zero(bmat):
            block_matrix[(i, j)] = bmat

    return block_matrix


def _calc_X0_inv(mat_X0, dim, n_iw):
    assert mat_X0.shape == (n_iw, dim, dim)

    mat_X0_inv = np.zeros((n_iw, dim, dim), dtype=complex)
    for iw in range(n_iw):
        mat_X0_inv[iw, :, :] = linalg.inv(mat_X0[iw])
    return mat_X0_inv


def _eigh_descending(mat):
    eigvals, eigvecs = linalg.eigh(mat)
    return eigvals[::-1], eigvecs[:, ::-1]


def _is_hermitian(mat):
    antihermite = (mat - mat.conjugate().T) / 2
    return _is_zero(antihermite)


def _hermitianize(mat):
    return (mat + mat.conjugate().T) / 2


def _load_chi_loc(h5, wb):
    chi_loc = h5.get(key=('chi_loc_in', wb))
    if chi_loc is False:
        print("  'chi_loc_in' not found. Compute 'chi_loc_in' from 'X_loc'.")
        # Compute chi_loc from X_loc
        beta = h5.get(key=('beta',))
        x_loc = h5.get(key=('X_loc', wb))
        if x_loc is False:
            sys.exit("Either chi_loc_in or X_loc is required.")
        chi_loc = {key: np.sum(array, axis=(2, 3))/beta for key, array in list(x_loc.items())}
    else:
        print("  'chi_loc_in' found")
    return chi_loc


def _copy_info(h5in, h5out):
    # read info
    keys = ['beta', 'block_name', 'inner_name']
    for key in keys:
        data = h5in.get(key)
        h5out.save(key=key, data=data)


# def calc_Iq_scl(prms : dict, h5_file : str, wb : np.ndarray, retain : int, remove : int, verbose=True):
def calc_Iq_scl(prms : dict):

    file_in = prms['input']
    file_out = prms['output']
    verbose = prms['verbose']

    print(f"\nLoad data from '{file_in!r}'")
    h5in = h5BSE(file_in)
    h5in.open('r')  # Keep HDF5 file open to improve performance. Close manually.

    beta = h5in.get(key=('beta'))
    print(f"  beta = {beta}")

    inner_name = h5in.get(key=('inner_name'))
    block_name = h5in.get(key=('block_name'))
    print(f"  inner_name: {inner_name}")
    print(f"  block_name: {block_name}")
    n_inner = len(inner_name)
    n_block = len(block_name)
    dim = n_inner * n_block

    # -------------------------------------------------------------
    # chi_loc

    print(f"\nLoad 'chi_loc'")
    wb = prms['wb']
    chi_loc = _load_chi_loc(h5in, wb)
    assert isinstance(chi_loc, dict)

    # Convert block matrix (dict) to ndarray
    mat_chi_loc = _convert_to_matrix(chi_loc, n_block, n_inner)
    assert mat_chi_loc.shape == (dim, dim)

    # hermitianize chi_loc
    if not _is_hermitian(mat_chi_loc):
        print("Warning: chi_loc is not hermitian.", file=sys.stderr)
    mat_chi_loc = _hermitianize(mat_chi_loc)

    # Get eigenvalues in descending order
    print("\nEigenvalues of chi_loc")
    eigs, u = _eigh_descending(mat_chi_loc)
    assert eigs.shape == (dim,)
    assert u.shape == (dim, dim)

    # TODO: retain, remove not necessary?

    # Number of eigenmodes retained
    # n_retain = retain
    # if n_retain is None:
    #     n_retain = dim - remove
    # print(f"  Number of eigenmodes retained = {n_retain}")

    for i, eig in enumerate(eigs):
        # if i == n_retain:
        #     print(f"  ------------------------- Remove eigenmodes below")
        print(f"  {i}: {eig}")

    # Dimensionality reduction
    # eigs = eigs[:n_retain]
    # u = u[:, :n_retain]
    # u_dag = u.conj().T
    # assert eigs.shape == (n_retain,)

    # -------------------------------------------------------------
    # X0_loc(iw)

    print(f"\nLoad 'X0_loc'")
    X0_loc = h5in.get(key=('X0_loc', wb))
    assert isinstance(X0_loc, dict)

    # Get n_iw (including negative freq)
    n_iw = X0_loc[list(X0_loc.keys())[0]].shape[2]  # (n_inner, n_inner, n_iw)
    print(f"  n_iw = {n_iw}")

    # TODO: iw_cutoff
    if prms['iw_cutoff'] != 0:
        raise NotImplementedError

    # Convert block matrix (dict) to ndarray
    mat_X0_loc = _convert_to_matrix_iw(X0_loc, n_block, n_inner, n_iw)
    assert mat_X0_loc.shape == (n_iw, dim, dim)

    # X0_loc^{-1}(iw)
    mat_X0_loc_inv = _calc_X0_inv(mat_X0_loc, dim, n_iw)
    assert mat_X0_loc_inv.shape == (n_iw, dim, dim)

    # -------------------------------------------------------------
    # phi(iw)

    wn = (2 * np.arange(-n_iw//2, n_iw//2) + 1) * np.pi / beta
    assert wn.shape == (n_iw,)
    assert np.isclose(wn[0], -wn[-1])  # wn_min == -wn_max
    print(f"  iw_max = {wn[-1]}")

    E_p = prms['E_plus']
    E_m = prms['E_minus']
    phi = 1 / (1j*wn + E_m) - 1 / (1j*wn - E_p)
    phi_sq = phi**2

    # -------------------------------------------------------------
    # X0_q(iw)

    # Get list of 'X0_q' data
    keylist = h5in.get_keylist_data(input_output='input')
    keylist_X0q = [key for key in keylist if key[0:2]==('X0_q', wb)]
    h5in.close()

    print(f"\n{len(keylist_X0q)} data of 'X0_q' found")
    print(f"Compute 'I_q_scl' and 'chi_q_scl' (existing data are overwritten)")

    print(f"\nResults are saved to {file_out!r}")
    h5out = h5BSE(file_out)
    h5out.open('a')  # Keep HDF5 file open to improve performance. Close manually.
    if file_in == file_out:
        h5in = h5out
    else:
        h5in.open('r')
        _copy_info(h5in, h5out)

    h5out.save(key=('chi_loc', wb), data=chi_loc)

    print(f"\nStart q-loop")
    for key in keylist_X0q:
        if verbose:
            print(f"Load data: key={key}")
        X0_q = h5in.get(key=key)
        assert isinstance(X0_q, dict)

        # Convert block matrix (dict) to ndarray
        mat_X0_q = _convert_to_matrix_iw(X0_q, n_block, n_inner, n_iw)
        assert mat_X0_q.shape == (n_iw, dim, dim)

        # Lambda_q(iw) = X0_loc^{-1}(iw) - X0_q^{-1}(iw)
        Lambda_q = mat_X0_loc_inv - _calc_X0_inv(mat_X0_q, dim, n_iw)
        assert Lambda_q.shape == (n_iw, dim, dim)

        # Compute I_q : sum over Matsubara freqs
        I_q = np.sum(Lambda_q * phi_sq[:, None, None], axis=0) / beta
        assert I_q.shape == (dim, dim)

        # Compute chi_q
        chi_q = mat_chi_loc @ linalg.inv(np.identity(dim) - I_q @ mat_chi_loc)

        # Convert ndarray to block matrix (dict)
        I_q_dict = _convert_to_dict(I_q, n_block, n_inner)
        chi_q_dict = _convert_to_dict(chi_q, n_block, n_inner)
        assert isinstance(I_q_dict, dict)
        assert isinstance(chi_q_dict, dict)

        _, _w, q = key
        assert _w == wb
        if verbose:
            print(f"Save data: key=('I_q_scl', {wb}, {q})")
            print(f"Save data: key=('chi_q_scl', {wb}, {q})")
        h5out.save(key=('I_q_scl', wb, q), data=I_q_dict)
        h5out.save(key=('chi_q_scl', wb, q), data=chi_q_dict)

    print("End q-loop")

    h5out.close()
    if file_in != file_out:
        h5in.close()


def read_params(file):
    print(f"\nRead file '{file}'")

    # check if file exist
    if not os.path.isfile(file):
        sys.exit(f"ERROR: File '{file}' not found.")

    config = configparser.ConfigParser()
    default_config = {
        'SCL': {
            'input' : 'dmft_bse.h5',
            'output' : 'dmft_bse.out.h5',
            'iw_cutoff' : 0,
            'wb' : 0,
            # 'remove' : 1,
            'verbose' : False,
        },
    }
    config.read_dict(default_config)
    config.read(file)

    prms = OrderedDict()
    prms['input'] = config.get('SCL', 'input')
    prms['output'] = config.get('SCL', 'output')
    prms['E_plus'] = config.getfloat('SCL', 'E_plus')
    prms['E_minus'] = config.getfloat('SCL', 'E_minus')
    prms['iw_cutoff'] = config.getint('SCL', 'iw_cutoff')
    prms['wb'] = config.getint('SCL', 'wb')
    # prms['remove'] = config.getint('SCL', 'remove')
    # prms['retain'] = config.getint('SCL', 'retain')
    prms['verbose'] = config.getboolean('SCL', 'verbose')

    # print(prms)
    for key, val in prms.items():
        print(f"  {key:10} = {val}")
    return prms


def main():
    P = argparse.ArgumentParser()
    P.add_argument('input_file', help="Input parameter file")
    P.add_argument('--verbose', action='store_true', help="verbose output")
    args = P.parse_args()
    print(args)

    prms = read_params(args.input_file)

    calc_Iq_scl(prms=prms)


if __name__ == "__main__":
    main()
