#!/usr/bin/env python3

import numpy as np
from fractions import Fraction
import sys
import configparser
import h5py
import os
import argparse
import ast


# Move q point into 0 <= q < 2*pi  (for integer index, 0 <= q < L)
def regularize(q, L):
    if q < 0:
        q += (-q // L + 1) * L
    return q % L


class GenQPath:
    """Generator of q-path file which is readable by BSE solver"""

    def __init__(self):
        self.L = None
        self.bvec = np.identity(3, dtype=float)
        self.qpoints = None

    def set_L_from_config(self, _filein_config):
        """Set L (system size) from config file"""
        print(f"\nRead k-grid from file {_filein_config!r}.")

        config = configparser.ConfigParser()
        config.read(_filein_config)

        # BSE format
        # if section [H0] exists, get L from config file
        if config.has_section('H0'):
            Lx = config.getint('H0', 'Lx')
            Ly = config.getint('H0', 'Ly')
            Lz = config.getint('H0', 'Lz')

            print("[H0]")
            print(" (Lx, Ly, Lz) = (%d, %d, %d)" %(Lx, Ly, Lz))
            self.L = [Lx,Ly,Lz]

        # DCore format
        elif config.has_section('model'):
            if config.has_option('model', 'nk'):
                nk0 = nk1 = nk2 = config.getint('model', 'nk')
            else:
                nk0 = config.getint('model', 'nk0')
                nk1 = config.getint('model', 'nk1')
                nk2 = config.getint('model', 'nk2')

            if config.has_option('model', 'lattice'):
                lattice = config.get('model', 'lattice')
                if lattice == 'square':
                    nk2 = 1

            print("[model]")
            print(" (nk0, nk1, nk2) = (%d, %d, %d)" % (nk0, nk1, nk2))
            self.L = [nk0, nk1, nk2]

        # try to find HDF5
        elif config.has_option('DMFT', 'dft_h5file'):
            dft_h5file = config.get('DMFT', 'dft_h5file')
            self.set_L_from_h5(dft_h5file)

        else:
            sys.exit(f"ERROR: File {_filein_config!r} does not contain information on k-grid.")

    def set_L_from_h5(self, _filein_h5):
        """set L (system size) from h5 file"""
        print(f"\nRead k-grid from file {_filein_h5!r}.")

        h5 = h5py.File(_filein_h5, 'r')
        self.L = h5['dft_input_chi/div'][()]
        assert isinstance(self.L, np.ndarray)
        print("L =", self.L)
        del h5
        assert len(self.L) == 3

    def set_bvec_from_config(self, _filein_config):
        """Set bvec (reciprocal vectors) from config file"""
        print(f"\nRead bvec (reciprocal vectors) from file {_filein_config!r}.")

        config = configparser.ConfigParser()
        config.read(_filein_config)

        # DCore format
        if config.has_option('model', 'bvec'):
            bvec_str = config.get('model', 'bvec')
            print("[model]")
            print("bvec =", bvec_str)
            bvec = ast.literal_eval(bvec_str)
            assert isinstance(bvec, list)

            self.bvec = np.array(bvec, dtype=float).T  # column vectors
            assert self.bvec.shape == (3,3)
            print("\nbvec (converted to column vectors) =\n", self.bvec)

        else:
            print("  bvec not found. Assume the identity matrix.")

    def set_qpoints(self, _filein):
        """Set qpoints (generator of q-path). This function should called after L is set"""
        print(f"\nRead q-points from {_filein!r} and convert to grid points.")

        def frac2i(s, _L):
            x = Fraction(s)*_L
            # print x, x.denominator
            if x.denominator != 1:
                sys.exit(f"ERROR: {s} is not on a grid point: {s}*{_L} = {x}.")
            return x.numerator

        self.qpoints = []
        with open(_filein, "r") as f:
            for line in f:
                if line[0] == '#':  # comment line
                    continue
                q = line.split()
                qint = np.array([frac2i(q[i], self.L[i]) for i in range(3)])
                qfloat = np.array([float(Fraction(q[i])) for i in range(3)])
                self.qpoints.append({'float': qfloat, 'int': qint, 'name': q[3]})

        # print qpoints
        print("q-path runs through (w.r.t. reciprocal primitive vectors)")
        for q in self.qpoints:
            print("", q['int'], q['name'])

    def gen_qpath(self, _fileout):
        """Generate q-path. This function should be called after all member variables are set"""

        # check if necessary quantities have been set
        for var in ['L', 'bvec', 'qpoints']:
            if self.__dict__[var] is None:
                raise ValueError("'%s' has not been set" % var)

        def qvec2str(q):
            # j = ( (q[0] * Ly) + q[1] ) * Lz + q[2]
            # return "%04d" %j
            q_reg = [regularize(q[i], self.L[i]) for i in range(3)]
            return "%02d.%02d.%02d" %(q_reg[0],q_reg[1],q_reg[2])

        def on_lattice(q_frac):
            for i in range(3):
                if q_frac[i].denominator != 1:
                    return False
            return True

        # output q-path
        with open(_fileout, 'w') as f:
            # the first point
            q0 = self.qpoints[0]
            x = 0
            q0['x'] = x
            print("%s %9.5f" % (qvec2str(q0['int']), x), q0['name'], file=f)

            # after the second point
            for q in self.qpoints[1:]:
                dif = q['int'] - q0['int']
                div = max( abs(dif[0]), abs(dif[1]), abs(dif[2]) )
                dq = np.array([Fraction(dif[0],div), Fraction(dif[1],div), Fraction(dif[2],div)])

                # compute distance in cartesian coordinate
                dq_float = q['float'] - q0['float']
                dq_cart = np.dot(self.bvec, dq_float)  # from primitive to cartesian
                dx = np.linalg.norm(dq_cart) / div

                qint = q0['int']
                for n in range(div-1):
                    qint = qint + dq
                    x += dx
                    if on_lattice(qint):
                        qvec = [qint[0].numerator, qint[1].numerator, qint[2].numerator]
                        print("%s %9.5f" %(qvec2str(qvec), x), file=f)

                x += dx
                q['x'] = x
                print("%s %9.5f" %(qvec2str(q['int']), x), q['name'], file=f)
                q0 = q

        print("\nFile '%s' generated" %_fileout)

        print("\nx-coordinate for plot")
        for q in self.qpoints:
            print("%9.5f %s" %(q['x'], q['name']))


def main():
    P = argparse.ArgumentParser()
    P.add_argument('file_param', help="config file or HDF5 file. Config file should have (Lx, Ly, Lz) in [H0] section, (nk0, nk1, nk2) or nk in [model] section (DCore interface), or dft_h5file in [DMFT] section.")
    P.add_argument('file_qpoints', help="generator of q-path")
    P.add_argument('-o', '--outfile', type=str, default='q_path.dat', help="Output file name")
    args = P.parse_args()
    # print(args)

    def file_check(_filename):
        if not os.path.exists(_filename):
            sys.exit(f"ERROR: File '{_filename}' not found.")
    file_check(args.file_param)
    file_check(args.file_qpoints)

    Q = GenQPath()
    if os.path.splitext(args.file_param)[1] == '.h5':  # check the extension
        Q.set_L_from_h5(args.file_param)  # read HDF5 file
    else:
        Q.set_L_from_config(args.file_param)  # read config file
        Q.set_bvec_from_config(args.file_param)

    Q.set_qpoints(args.file_qpoints)
    Q.gen_qpath(args.outfile)


if __name__ == '__main__':
    main()
