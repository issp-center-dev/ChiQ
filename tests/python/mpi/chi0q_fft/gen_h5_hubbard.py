

from hk_converter_chi import * # In the same directory
from itertools import *
import numpy as np
import os
import sys

#=============================================================================
# set input file

argvs = sys.argv
argc = len(argvs)
if argc != 2:
    raise Exception("How to use: python3 %s file_param" %argvs[0])

file_param=argvs[1]
print("file_param =", file_param)
if not os.path.exists(file_param):
    raise Exception("File '%s' not exists" %file_param)

#=============================================================================
# set paramters

import configparser

default_config = {
    't' : '1',
    'dft_h5file' : 'hubbard.h5',
    'n_orb' : '1',
    'l_orb' : '0',
    'occup' : '1.0',
    'n_sublattice' : '1',
    'stag_pot' : '0.0',
}
config = configparser.ConfigParser(default_config)
config.read(file_param)

# read parameters
block = 'H0'
print("[%s]" % block)
params = {}
for key in ['Lx', 'Ly', 'Lz', 't', 'n_orb', 'l_orb', 'n_sublattice', 'stag_pot']:
    params[key] = config.get(block, key)
    print("", key, "=", params[key])
# print(params)
Lx = int(params['Lx'])
Ly = int(params['Ly'])
Lz = int(params['Lz'])
t = float(params['t'])
n_orb = int(params['n_orb'])
l_orb = int(params['l_orb'])
n_sublattice = int(params['n_sublattice'])
stag_pot = float(params['stag_pot'])
del params

# read parameters
block = 'DMFT'
print("[%s]" % block)
params = {}
for key in ['dft_h5file', 'occup']:
    params[key] = config.get(block, key)
    print("", key, "=", params[key])
dft_h5file = params['dft_h5file']
occup = float(params['occup'])
del params

# check parameters
def is_even_or_one(_n):
    return _n % 2 == 0 or _n==1
if n_sublattice == 1:
    assert stag_pot == 0.
elif n_sublattice == 2:
    # for AFM
    assert is_even_or_one(Lx)
    assert is_even_or_one(Ly)
    assert is_even_or_one(Lz)
else:
    raise ValueError("not implemented")

#=============================================================================

file_name, ext = os.path.splitext(dft_h5file)

def get_cosk(L):
    if L==1:
        return [0,]
    return [ np.cos(2.0*i*np.pi/L) for i in range(L) ]

coskx = get_cosk(Lx)
cosky = get_cosk(Ly)
coskz = get_cosk(Lz)

def get_ek(_i, _j, _k):
    # e(k)
    return -2.0 * t * (coskx[_i] + cosky[_j] + coskz[_k])

def get_ekq(_i, _j, _k):
    # e(k+Q)
    i2 = (_i + Lx // 2) % Lx
    j2 = (_j + Ly // 2) % Ly
    k2 = (_k + Lz // 2) % Lz
    return get_ek(i2, j2, k2)

def in_BZ(_i, _j, _k):
    """return True is (i,j,k) is in the BZ"""
    if _i < Lx // n_sublattice:
        return True
    else:
        return False

with open(file_name, 'w') as f:
    # N = Lx*Ly*Lz/n_sublattice
    N = Lx*Ly*Lz  # sum up all k-points in the original (unfolded) BZ
    print(N, file=f)   # number of k-points
    print(occup*n_sublattice, file=f)      # Electron density
    print(n_sublattice, file=f)            # number of total atomic shells
    for a in range(n_sublattice):
        print("%d %d %d %d" %(a, a, l_orb, n_orb), file=f)  # iatom, isort, l, dimension
    print(n_sublattice, file=f)            # number of correlated shells
    for a in range(n_sublattice):
        print("%d %d %d %d 0 0" % (a, a, l_orb, n_orb), file=f)  # iatom, isort, l, dimension, SO, irep
    for a in range(n_sublattice):
        # print("1 1", file=f)               # Number of ireps, dimension of irep
        print("1 %d" %n_orb, file=f)               # Number of ireps, dimension of irep
        # TODO: check results for multiorbital sample

    count=0
    for i,j,k in product(list(range(Lx)),list(range(Ly)),list(range(Lz))):
        # if not in_BZ(i, j, k):
        #     continue

        ek = get_ek(i, j, k)  # e(k)
        ekq = get_ekq(i, j, k)  # e(k+Q)

        if n_sublattice==1:
            hmat = np.array([[ek,],])
        elif n_sublattice==2:
            # Hamiltonian matrix for 2-sublattice
            hmat = np.array([[(ek + ekq) / 2 + stag_pot / 2, (ek - ekq) / 2],
                             [(ek - ekq) / 2, (ek + ekq) / 2 - stag_pot / 2]])

        for func in [np.real, np.imag]:
            for a1, a2 in product(list(range(n_sublattice)), repeat=2):
                for o1,o2 in product(list(range(n_orb)), repeat=2):
                    print("%.20e" %func(hmat[a1,a2]) if o1==o2 else 0.0, file=f)
        count += 1
    assert count == N

Converter = HkConverterChi(filename = file_name)
Converter.convert_dft_input()

os.remove(file_name)

# for chi0 calculations
div = [Lx, Ly, Lz]  # all k-points in the original (unfolded) BZ
# kpoints = []
# k2i = np.empty(div, dtype=int)
# for m, [i,j,k] in enumerate(product(range(Lx),range(Ly),range(Lz))):
#     kpoints.append([i,j,k])
#     k2i[i,j,k] = m

# Converter.convert_dft_input_chi(div, kpoints, k2i)
Converter.convert_dft_input_chi(div)
