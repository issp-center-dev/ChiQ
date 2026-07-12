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

    shutil.copy('ref/dmft_bse.h5', './')  # overwrite if exists

    # run
    assert os.system("chiq_main.py bse.toml") == 0
    assert os.system("chiq_post.py bse.toml") == 0

    # compare results (take rounding errors into account)
    chi0q = read_chiq('chi0_q_eigen.dat')
    chi0q_ref = read_chiq('ref/chi0_q_eigen.dat')
    assert np.allclose(chi0q, chi0q_ref)

    chiq = read_chiq('chi_q_eigen.dat')
    chiq_ref = read_chiq('ref/chi_q_eigen.dat')
    assert np.allclose(chiq, chiq_ref)


def test_chiq_numpy_backend(request):
    org_dir = os.getcwd()
    os.chdir(request.fspath.dirname)

    try:
        shutil.copy('ref/dmft_bse.h5', './')  # overwrite if exists

        # generate a numpy-backend variant of bse.toml
        with open('bse.toml') as f:
            toml_text = f.read()
        toml_text_numpy = toml_text.replace(
            '[chiq_main]', '[chiq_main]\nbackend = "numpy"'
        )
        assert toml_text_numpy != toml_text
        with open('bse_numpy.toml', 'w') as f:
            f.write(toml_text_numpy)

        # run
        assert os.system("chiq_main.py bse_numpy.toml") == 0
        assert os.system("chiq_post.py bse_numpy.toml") == 0

        # compare results against the same reference the cpp backend matches
        chi0q = read_chiq('chi0_q_eigen.dat')
        chi0q_ref = read_chiq('ref/chi0_q_eigen.dat')
        assert np.allclose(chi0q, chi0q_ref)

        chiq = read_chiq('chi_q_eigen.dat')
        chiq_ref = read_chiq('ref/chi_q_eigen.dat')
        assert np.allclose(chiq, chiq_ref)
    finally:
        os.chdir(org_dir)
