#!/usr/bin/env python3


import numpy as np
import os
from itertools import product
import argparse
import configparser

from chiq.h5bse import h5BSE


def print_dtype_shape(np_array):
    assert isinstance(np_array, np.ndarray)
    print("  dtype =", np_array.dtype, "  shape =", np_array.shape)


def main():
    P = argparse.ArgumentParser()
    P.add_argument('file_param', help="config file")
    args = P.parse_args()
    # print args

    def file_check(_filename):
        if not os.path.exists(_filename):
            raise Exception("File '%s' not exists" %_filename)
    file_check(args.file_param)

    # -------------------------------------------------------------------
    # set parameters

    default_config = {
        'ext' : 'txt',
        'h5_file' : 'dmft_bse.h5',
        'groupname' : '',
        'qpoints_saved' : 'quadrant',
    }
    config = configparser.ConfigParser(default_config)
    config.read(args.file_param)

    # [H0]
    Lx = int(config.get('H0', 'Lx'))
    Ly = int(config.get('H0', 'Ly'))
    Lz = int(config.get('H0', 'Lz'))

    # [DMFT]
    beta = float(config.get('DMFT', 'beta'))

    # [TXT2BSE]
    basename = config.get('TXT2BSE', 'basename')
    ext = config.get('TXT2BSE', 'ext')
    h5_file = config.get('TXT2BSE', 'h5_file')
    groupname = config.get('TXT2BSE', 'groupname')
    qpoints_saved = config.get('TXT2BSE', 'qpoints_saved')

    # add '.' if not included
    if ext[0] != '.':
        ext = '.' + ext
    # print ext

    # -------------------------------------------------------------------
    # save to hdf5
    BS = h5BSE(h5_file, groupname)

    # info
    block_name = ['sp']
    inner_name = ['0-0']
    print("beta =", beta)
    print("block_name =", block_name)
    print("inner_name =", inner_name)
    BS.save( key=('beta'), data=beta )
    BS.save( key=('block_name'), data=block_name )
    BS.save( key=('inner_name'), data=inner_name )


    def convert_txt_to_h5(key, postfix, flag_complex):
        print("\n%s" %key)
        filename = basename + postfix + ext
        if os.path.exists(filename):
            print(" reading file '%s'" %filename)
            if flag_complex:
                data = np.loadtxt(filename).view(complex)
            else:
                data = np.loadtxt(filename)
            # print data
            data_reshape = np.array([[data]], dtype=np.complex128)
            print_dtype_shape(data_reshape)
            BS.save( key=(key,0), data={ (0,0) : data_reshape } )
        else:
            print(" file '%s' not found" %filename)

    # chi0_loc
    convert_txt_to_h5(key='chi0_loc', postfix='-chi0loc', flag_complex=False)

    # chi_loc
    convert_txt_to_h5(key='chi_loc', postfix='-chiloc', flag_complex=False)

    # X0_loc(iw)
    convert_txt_to_h5(key='X0_loc', postfix='-x0loc', flag_complex=True)

    # X_loc(iw,iw')
    convert_txt_to_h5(key='X_loc', postfix='-xloc', flag_complex=True)


    # X0_q(iw)
    print("\nX0_q(iw)")

    q_int2str = { i : '%02d.%02d.%02d' %(qx,qy,qz) for i,(qx,qy,qz) in enumerate(product(list(range(Lx)),list(range(Ly)),list(range(Lz)))) }
    # print q_int2str

    q_save = [False]*(Lx*Ly*Lz)
    print(" qpoints_saved = %s" % qpoints_saved)
    if qpoints_saved == 'quadrant':  # q in the first quadrant
        for kx,ky,kz in product(list(range(Lx/2+1)), list(range(Ly/2+1)), list(range(Lz/2+1))):
            q = (kx * Ly + ky ) * Lz + kz
            q_save[q] = True
    elif qpoints_saved == 'fbz':
        q_save = [True] * (Lx * Ly * Lz)
    else:
        print("ERROR: qpoints_saved=%s is not supported" % qpoints_saved)
        exit(1)

    # print q_save

    filename = basename + '-x0q' + ext
    if os.path.exists(filename):
        data = np.loadtxt(filename, dtype=str)
        # print data
        # print type(data)
        print_dtype_shape(data)
        count = 0
        for line in data:
            # print line
            q = int(line[0])
            # print q
            if q_save[q]:
                data = line[1:].astype(np.float64).view(complex)
                X0_q = np.array([[data]])
                # print data
                BS.save( key=('X0_q', 0, q_int2str[q]), data={ (0,0) : X0_q } )
                count += 1
        print(" %d data saved" % count)
    else:
        print(" file '%s' not found" %filename)


if __name__ == '__main__':
    main()
