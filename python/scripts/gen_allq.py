#!/usr/bin/env python3

import os
import sys
import argparse
from itertools import product

from gen_qpath import GenQPath # in the same directory


class GenQfbz(GenQPath):

    def gen_fbz(self, _fileout):
        """Generate all q in FBZ"""

        def qvec2str(q):
            return "%02d.%02d.%02d" %(q[0],q[1],q[2])

        # output all q
        with open(_fileout, 'w') as f:
            for w in self.wlist:
                for qx, qy, qz in product(list(range(self.L[0])), list(range(self.L[1])), list(range(self.L[2]))):
                    qvec = [qx, qy, qz]
                    print("%s %s" % (w, qvec2str(qvec)), file=f)

        print("\nFile '%s' generated." %_fileout)


def main():
    P = argparse.ArgumentParser()
    P.add_argument('file_param', help="config file or HDF5 file. Config file should have either (Lx, Ly, Lz) in [H0] section or dft_h5file in [DMFT] section.")
    P.add_argument('-o', '--outfile', type=str, default='q_fbz.dat', help="Output file name")
    args = P.parse_args()
    # print(args)

    def file_check(_filename):
        if not os.path.exists(_filename):
            sys.exit(f"ERROR: File '{_filename}' not found.")
    file_check(args.file_param)

    Q = GenQfbz()
    if os.path.splitext(args.file_param)[1] == '.h5':  # check the extension
        Q.set_L_from_h5(args.file_param)  # read HDF5 file
    else:
        Q.set_L_from_config(args.file_param)  # read config file

    Q.set_wlist([0,])
    Q.gen_fbz(args.outfile)


if __name__ == '__main__':
    main()
