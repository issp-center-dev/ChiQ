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

    shutil.copy('ref/dmft_bse.out.h5', './')  # overwrite if exists

    # run
    assert os.system("chiq_post.py bse.toml") == 0

    # compare results (take rounding errors into account)
    chi0q = read_chiq('chi0_q_linear_combination.dat')
    chi0q_ref = read_chiq('ref/chi0_q_linear_combination.dat')
    assert np.allclose(chi0q, chi0q_ref)

    chiq = read_chiq('chi_q_linear_combination.dat')
    chiq_ref = read_chiq('ref/chi_q_linear_combination.dat')
    assert np.allclose(chiq, chiq_ref)
