import os
import shutil


def test_gk_mode_reproduces_saved_x0(request):
    os.chdir(request.fspath.dirname)
    for f in ('hubbard.h5', 'x0_save.h5', 'x0_gk.h5', 'x0_gk_temp.h5', 'G_k_iw_temp.h5'):
        if os.path.exists(f):
            os.remove(f)
    for f in ('gen_h5_hubbard.py', 'hk_converter_chi.py'):
        shutil.copy(os.path.join('..', 'chi0q_fft', f), f)

    assert os.system("python3 gen_h5_hubbard.py param.in") == 0  # generate hubbard.h5
    # SumkDFT imports mpi4py, so it should be executed with mpirun
    # allow 2 ranks on small CI runners (OpenMPI; ignored by other MPI implementations)
    os.environ.setdefault("OMPI_MCA_rmaps_base_oversubscribe", "1")
    assert os.system("mpirun -np 2 python3 x0_gk_modes.py param.in") == 0
    assert os.system("python3 compare_modes.py") == 0
