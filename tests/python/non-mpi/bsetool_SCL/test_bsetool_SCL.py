import os
import shutil
import numpy as np


def read_chiq(filename):
    # drop the first two columns
    array = np.loadtxt(filename, dtype=str)[:,2:].astype('float')
    return array

def test_chiq(request):
    org_dir = os.getcwd()
    os.chdir(request.fspath.dirname)

    shutil.copy('../bsetool_BSE/ref/dmft_bse.h5', './')  # overwrite if exists

    # run
    assert os.system("calc_Iq_scl.py scl_2pole.in") == 0
    assert os.system("bse_post.py bse.toml") == 0

    # compare results (take rounding errors into account)

    chiq = read_chiq('chi_q_scl_eigen.dat')
    chiq_ref = read_chiq('ref/chi_q_scl_eigen.dat')
    assert np.allclose(chiq, chiq_ref)
