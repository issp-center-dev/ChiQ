#!/usr/bin/env python3

import os
import argparse
from itertools import product
import copy
import numpy as np
from scipy import fftpack

from chiq.h5bse import h5BSE
from gen_qpath import GenQPath # in the same directory


def isnonzero(array):
    """
    check if np.ndarray is nonzero
    """
    return np.linalg.norm(array) > 1e-10


class BSEFFT(object):

    def __init__(self, nk, file_bse):
        assert len(nk) == 3
        self._nk = copy.copy(nk)

        # set in get_data()
        self.data = None
        self.h5bse = None
        self.n_block = None
        self.n_inner = None

        # set in reform_data()
        self.data_q = None

        # set in fft_data()
        self.data_r = None

        self.h5bse = h5BSE(file_bse)

    def get_data(self, input_dname='I_q', input_w=0):
        print("\nReading data...")

        self.h5bse.open('r')  # Keep HDF5 file open to improve performance. Close manually.

        self.data = {}
        for qx, qy, qz in product(list(range(self._nk[0])), list(range(self._nk[1])), list(range(self._nk[2]))):
            str_q = "%02d.%02d.%02d" % (qx, qy, qz)
            key = (input_dname, input_w, str_q)
            data = self.h5bse.get(key)
            self.data[(qx, qy, qz)] = data

        print("  done len=%d" % len(self.data))

        block_names = self.h5bse.get("block_name")
        inner_names = self.h5bse.get("inner_name")
        print("  block_names =", block_names)
        print("  inner_names =", inner_names)
        self.n_block = len(block_names)
        self.n_inner = len(inner_names)

        self.h5bse.close()

    def reform_data(self):
        print("\nReforming data...")

        # [qx, qy, qz, i, j]
        #     i, j for spin-orbital
        qvec = (self._nk[0], self._nk[1], self._nk[2])
        mat_size = self.n_block * self.n_inner
        shape_full = qvec + (mat_size, mat_size)
        self.data_q = np.zeros(shape_full, dtype=complex)
        print(" ", shape_full)

        for qx, qy, qz in product(list(range(self._nk[0])), list(range(self._nk[1])), list(range(self._nk[2]))):
            data = self.data[(qx, qy, qz)]

            for block, mat in list(data.items()):
                # position to be inserted
                b1, b2 = block
                i = b1 * self.n_inner
                j = b2 * self.n_inner
                self.data_q[qx, qy, qz, i:i+self.n_inner, j:j+self.n_inner] = mat[:, :]

        self.data = None  # release memory

    def fft_data(self):
        print("\nPerforming FFT...")

        # [qx, qy, qz, i, j] -> [i, j, qx, qy, qz]
        temp_q = self.data_q.transpose((3, 4, 0, 1, 2))
        temp_r = np.zeros(temp_q.shape, dtype=complex)

        mat_size = temp_q.shape[0]
        for i, j in product(list(range(mat_size)), repeat=2):
            # FFT: q -> r
            temp_r[i, j] = fftpack.ifftn(temp_q[i, j], overwrite_x=True)

        # [i, j, rx, ry, rz] -> [rx, ry, rz, i, j]
        self.data_r = temp_r.transpose((2, 3, 4, 0, 1))

        self.data_q = None  # release memory

    def save_data(self, output_dname='I_r', output_w=0):
        print("\nSaving data into HDF5 file...")

        self.h5bse.open('a')  # Keep HDF5 file open to improve performance. Close manually.

        for rx, ry, rz in product(list(range(self._nk[0])), list(range(self._nk[1])), list(range(self._nk[2]))):
            data = {}
            for b1, b2 in product(list(range(self.n_block)), repeat=2):
                i = b1 * self.n_inner
                j = b2 * self.n_inner
                mat = self.data_r[rx, ry, rz, i:i+self.n_inner, j:j+self.n_inner]
                if isnonzero(mat):
                    data[(b1, b2)] = mat

            # If data is empty, set zero matrix at (0,0) block to avoid error in bse_tool.
            if not data:
                data[(0, 0)] = np.zeros((self.n_inner, self.n_inner), dtype=complex)

            str_r = "%02d.%02d.%02d" % (rx, ry, rz)
            key = (output_dname, output_w, str_r)

            self.h5bse.save(key, data)

        self.h5bse.close()


def main():
    P = argparse.ArgumentParser()
    P.add_argument('file_param', help="config file or HDF5 file. Config file should have (Lx, Ly, Lz) in [H0] section, (nk0, nk1, nk2) or nk in [model] section (DCore interface), or dft_h5file in [DMFT] section.")
    P.add_argument('-f', '--file_bse', default='dmft_bse.out.h5', help="input/output HDF5 file containing BSE data (default: 'dmft_bse.out.h5')")
    P.add_argument('--input_dname', default='I_q', help="input data name in HDF5 file (default: 'I_q')")
    P.add_argument('--output_dname', default='I_r', help="output data name in HDF5 file (default: 'I_r')")
    P.add_argument('-w', type=int, default=0, help="input/output frequency (default: 0)")

    args = P.parse_args()
    print(args)

    h5file = args.file_bse

    def file_check(_filename):
        if not os.path.exists(_filename):
            raise Exception("File '%s' not exists" % _filename)
    file_check(args.file_param)

    Q = GenQPath()
    if os.path.splitext(args.file_param)[1] == '.h5':  # check the extension
        Q.set_L_from_h5(args.file_param)  # read HDF5 file
    else:
        Q.set_L_from_config(args.file_param)  # read config file

    bse = BSEFFT(Q.L, h5file)
    bse.get_data(args.input_dname, args.w)
    bse.reform_data()
    bse.fft_data()
    bse.save_data(args.output_dname, args.w)


if __name__ == '__main__':
    main()
