#!/usr/bin/env python3

import numpy as np
from itertools import product
import os
import math
import argparse

# from bse_tools.h5bse import h5BSE
from chiq.h5bse import h5BSE

if __name__ == '__main__':

    print("\n-----START")
    print(__file__)

    P = argparse.ArgumentParser()
    P.add_argument('-f', default='dmft_bse.h5', help="h5 filename")
    P.add_argument('--scl', choices=['svd', '1dhalf', '2pole'], required=True)
    args = P.parse_args()
    print(args)

    h5_file = args.f

    groupname_list = { 'svd' : 'scl_svd', '1dhalf' : 'scl_1dhalf', '2pole' : 'scl_2pole' }
    groupname = groupname_list[args.scl]

    # get Phi_L kept in a subgroup
    BS_exact = h5BSE(h5_file, groupname=groupname)
    phi_dict = BS_exact.get( key=('Phi', 0) )

    # overwrite Phi_L in DB
    BS = h5BSE(h5_file)
    BS.save( key=('Phi', 0), data=phi_dict )

    print("Phi data was overwritten with those in group '%s'" %groupname)
    print("-----END")
