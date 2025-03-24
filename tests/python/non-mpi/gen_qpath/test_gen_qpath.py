import os
import filecmp

def test_gen_qpath(request):
    org_dir = os.getcwd()
    os.chdir(request.fspath.dirname)

    assert os.system("gen_qpath.py param.in qpath.in") == 0
    assert filecmp.cmp('q_path.dat', 'ref/q_path.dat')
