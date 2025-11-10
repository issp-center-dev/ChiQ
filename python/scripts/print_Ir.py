#!/usr/bin/env python3

import numpy as np
from scipy import linalg
import sys
import os
import argparse
import configparser
from itertools import product
from collections import OrderedDict
import matplotlib.pyplot as plt

from chiq.h5bse import h5BSE
from chiq.matrix_dict import MatrixDict


def _is_zero(mat, tol=1e-8):
    return np.all(np.absolute(mat) < tol)

def _is_real(mat, tol=1e-8):
    return _is_zero(np.imag(mat), tol)

def _is_hermitian(mat, tol=1e-8):
    antihermite = (mat - mat.conjugate().T) / 2
    return _is_zero(antihermite, tol)

def _hermitianize(mat):
    return (mat + mat.conjugate().T) / 2


def _convert_to_matrix(block_matrix, n_block, n_inner):
    assert isinstance(block_matrix, dict)
    dim = n_inner * n_block
    mat = np.zeros((n_block, n_inner, n_block, n_inner), dtype=complex)
    for block, bmat in block_matrix.items():
        i, j = block
        mat[i, :, j, :] = bmat
    return mat.reshape((dim, dim))


def convert_Rstr_to_Rvec(Rstr):
    Rvec = np.fromstring(Rstr, sep='.', dtype=int)  # 00.01.02 -> [0,1,2]
    return Rvec


def convert_frac_to_cart(R, frac_coord, avecs):
    assert R.shape == (3,)  # lattice vector in fractional coordinate
    assert frac_coord.shape == (3,)  # atomic position in fractional coordinate
    assert avecs.shape == (3, 3)

    #  <r|a> <a|x>  (<r|a> is real)
    return (R + frac_coord) @ avecs


def print_matrix(mat, tol=1e-12, **args):
    def chop(x, tol=tol):
        return f" 0           " if abs(x) < tol else f"{x:>13.6e}"
    assert mat.ndim == 2
    ncol, nrow = mat.shape
    for i in range(ncol):
        row_str = ""
        for j in range(nrow):
            row_str += f"{chop(mat[i,j].real)} {chop(mat[i,j].imag)} "
        print(row_str, **args)


def plot(x, y, basename, xmax):
    fig, ax = plt.subplots()
    ax.plot(x, y, '.')
    ax.set_xlabel(r"$|r|$")
    ax.set_ylabel(r"$I(r)$")
    ax.set_ylim(None, None)
    ax.set_axisbelow(True)
    ax.axhline(y=0, color="k", zorder=-1)  # show yzero

    # all range
    figname1 = basename + "_all.pdf"
    ax.set_xlim(0, None)
    fig.savefig(figname1)
    print(f"  '{figname1}'")

    # zoomed
    figname2 = basename + "_zoom.pdf"
    ax.set_xlim(0, xmax)
    fig.savefig(figname2)
    print(f"  '{figname2}'")

    # save data
    dat_name = basename + ".dat"
    np.savetxt(dat_name, np.column_stack((x, y)))
    print(f"  '{dat_name}'")


def print_Ir(prms : dict):

    file_in = prms['input']
    verbose = prms['verbose']
    n_sites = prms['n_sites']
    n_bases = prms['n_bases']
    sizes = prms['sizes']
    chop = prms['chop']

    output_dir = prms['output_dir']
    os.makedirs(output_dir, exist_ok=True)

    # -------------------------------------------------------------
    # Load data from HDF5 file

    print(f"\nLoad data from '{file_in}'")
    h5in = h5BSE(file_in)
    h5in.open('r')  # Keep HDF5 file open to improve performance. Close manually.

    inner_name = h5in.get(key=('inner_name'))
    block_name = h5in.get(key=('block_name'))
    print(f"  inner_name: {inner_name}")
    print(f"  block_name: {block_name}")
    n_inner = len(inner_name)
    n_block = len(block_name)
    dim = n_inner * n_block

    # -------------------------------------------------------------
    # Load data from text files

    # Load data from files
    print(f"\nLoad data from '{prms['file_avecs']}'")
    avecs = np.loadtxt(prms['file_avecs'])
    print(avecs.shape, avecs.dtype)
    assert avecs.shape == (3, 3), f"avecs.shape={avecs.shape} must be (3, 3)"

    print(f"\nLoad data from '{prms['file_coords']}'")
    coords = np.loadtxt(prms['file_coords'])
    print(coords.shape, coords.dtype)
    assert coords.shape == (n_sites, 3), f"coords.shape={coords.shape} must be ({n_sites}, 3)"

    print(f"\nLoad data from '{prms['file_bases']}'")
    bases = np.loadtxt(prms['file_bases'])
    print(bases.shape, bases.dtype)
    bases = bases.view(np.complex128)
    print(bases.shape, bases.dtype)
    assert bases.shape == (n_bases, dim), f"bases.shape={bases.shape} must be ({n_bases}, {dim})"

    # U = <m|a>
    U = bases.T
    assert U.shape == (dim, n_bases)

    # -------------------------------------------------------------
    # Site

    # sites
    sites = np.fromstring(prms['sites'], dtype=int, sep=' ')
    assert sites.shape == (n_bases,), f"sites.shape={sites.shape} must be ({n_bases},)"

    print(f"\nParse sites")
    sites_set = set(sites)
    # print(f"  {sites_set}")
    assert sites_set == set(range(n_sites)), f"sites must composed of 0, 1, ..., {n_sites-1}."
    assert len(sites_set) == n_sites, f"n_sites={n_sites}"

    # mapping from site to basis index
    site2idx = []
    for s in range(n_sites):
        site2idx.append(np.array([i for i, site in enumerate(sites) if site==s], dtype=int))

    print(f" site: basis indices")
    for s in range(n_sites):
        print(f" {s:4}: {site2idx[s]}")

    # -------------------------------------------------------------
    # Rvecs (lattice vectors in fractional coordinate)

    wb = prms['wb']
    type = 'I_r_scl'

    # Get list of 'I_r' data
    keylist = h5in.get_keylist_data(input_output='output')
    # print(keylist)
    keylist_Ir = [key for key in keylist if key[0:2]==(type, wb)]
    print(f"\n{len(keylist_Ir)} data of '{type}' found")

    # Make R vectors
    Rvecs = np.array([convert_Rstr_to_Rvec(key[2]) for key in keylist_Ir])

    # Get system size
    if sizes == 'auto':
        sizes = Rvecs.max(axis=0) + 1
        print(f" System size = {sizes}")
        assert sizes.shape == (3,)

    # Shift Rvecs into the first Brillouin zone
    Rvecs = np.where(Rvecs > sizes[None, :]/2, Rvecs - sizes[None, :], Rvecs)

    # -------------------------------------------------------------
    # R-loop

    print(f"\nResults are saved to 'I_r_site*.dat' (*=0~{n_sites-1})")
    fs = []
    for site in range(n_sites):
        filename = f"I_r_site{site}.dat"
        fs.append(open(filename, 'w'))
        print(f"# I(R=0, r_{site}; R, r_site)", file=fs[-1])

    print(f"\nStart R-loop")
    R0 = np.zeros(3, dtype=int)
    dists = []
    data4plot = []
    for key, Rvec in zip(keylist_Ir, Rvecs):
        # print(key, Rvec)
        if verbose:
            print(f"Load data: key={key}")

        I_r = h5in.get(key=key)
        assert isinstance(I_r, dict)

        # Convert block matrix (dict) to ndarray
        mat_I_r = _convert_to_matrix(I_r, n_block, n_inner)
        assert mat_I_r.shape == (dim, dim)

        # Convert basis
        mat_I_r = np.einsum('am,mn,nb->ab', U.conj().T, mat_I_r, U)
        assert mat_I_r.shape == (n_bases, n_bases)

        # Decompose into sites
        for s1, s2 in product(range(n_sites), repeat=2):
            mat_sh = mat_I_r[np.ix_(site2idx[s1], site2idx[s2])]

            r1 = coords[s1]
            r2 = coords[s2]
            r_diff = convert_frac_to_cart(Rvec, r2, avecs) - convert_frac_to_cart(R0, r1, avecs)
            dist = np.linalg.norm(r_diff)

            print(f"# R={Rvec} site={s2} r_diff={r_diff} dist={dist:.4f}", file=fs[s1])
            print_matrix(mat_sh, tol=chop, file=fs[s1])

            dists.append(dist)
            data4plot.append(mat_sh.reshape(-1))

    print("End R-loop")

    for f in fs:
        f.close()
    h5in.close()

    # Convert to ndarray
    dists = np.array(dists)
    data4plot = np.array(data4plot)
    # print(dists.shape)
    # print(data4plot.shape)
    assert dists.shape[0] == data4plot.shape[0]

    # Analyze distances
    dists_unique = np.unique(np.round(dists, decimals=10))  # sorted in lexicographic order
    print("\nDistances:")
    print(dists_unique[:10])
    print(f"  min: {dists_unique[1]:.6e}")
    print(f"  max: {dists_unique[-1]:.6e}")
    print(f"  number of different distances: {len(dists_unique)}:")

    # Analyze I(r)
    print("\nI(r):")
    print(f"  is real: {_is_real(data4plot)}")
    max_imag = np.max(np.abs(np.imag(data4plot)))
    print(f"  max imag: {max_imag:.4e}")

    data_unique = np.unique(np.round(np.real(data4plot), decimals=10))  # sorted in lexicographic order
    print("\nLargest antiferroic (<0) values:")
    print(data_unique[:10])  # print first 12 values
    print("\nLargest ferroic (>0) values:")
    print(data_unique[:-10:-1])  # print last 12 values

    # Plot
    print("\nPlot")
    # up to 9-th neighbor
    plot(x=dists, y=np.real(data4plot), basename=os.path.join(output_dir, "Ir_dist"), xmax=dists_unique[10])


def read_params(file):
    print(f"\nRead file '{file}'")

    # check if file exist
    if not os.path.isfile(file):
        sys.exit(f"ERROR: File '{file}' not found.")

    config = configparser.ConfigParser()
    default_config = {
        'Ir': {
            'input' : 'dmft_bse.out.h5',
            'output_dir' : './',
            'wb' : 0,
            'verbose' : False,
            'sizes' : 'auto',
            'chop' : 1e-12,
        },
    }
    config.read_dict(default_config)
    config.read(file)

    prms = OrderedDict()
    # optional
    prms['input'] = config.get('Ir', 'input')
    prms['output_dir'] = config.get('Ir', 'output_dir')
    prms['wb'] = config.getint('Ir', 'wb')
    prms['verbose'] = config.getboolean('Ir', 'verbose')
    prms['sizes'] = config.get('Ir', 'sizes')
    prms['chop'] = config.getfloat('Ir', 'chop')

    # mandatory
    prms['type'] = config.get('Ir', 'type')
    prms['n_bases'] = config.getint('Ir', 'n_bases')
    prms['n_sites'] = config.getint('Ir', 'n_sites')
    prms['sites'] = config.get('Ir', 'sites')
    prms['file_avecs'] = config.get('Ir', 'file_avecs')
    prms['file_coords'] = config.get('Ir', 'file_coords')
    prms['file_bases'] = config.get('Ir', 'file_bases')

    # print(prms)
    for key, val in prms.items():
        print(f"  {key:11} = {val}")
    return prms


def main():
    P = argparse.ArgumentParser()
    P.add_argument('input_file', help="Input parameter file")
    P.add_argument('--verbose', action='store_true', help="verbose output")
    args = P.parse_args()
    print(args)

    prms = read_params(args.input_file)

    print_Ir(prms=prms)


if __name__ == "__main__":
    main()
