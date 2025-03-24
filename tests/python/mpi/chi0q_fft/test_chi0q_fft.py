import os
import numpy as np


def remove_if_exists(filename):
    if os.path.exists(filename):
        os.remove(filename)


def test_compare_FFT_and_sum(request):
    org_dir = os.getcwd()
    os.chdir(request.fspath.dirname)

    remove_if_exists('hubbard.h5')
    remove_if_exists('chi0q_fft.h5')
    remove_if_exists('chi0q_sum.h5')

    # run
    assert os.system("python3 gen_h5_hubbard.py param.in") == 0  # generate hubbard.h5
    assert os.system("gen_qpath.py param.in qpath.in") == 0  # generate q_path.dat

    # SumkDFT imports mpi4py, so it should be executed with mpirun
    assert os.system("mpirun -np 1 python3 chi0q_fft_sum.py param.in") == 0

    # compare
    assert os.system("h5diff --delta=1e-12 chi0q_fft.h5 chi0q_sum.h5") == 0
