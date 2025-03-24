import os
import numpy
from chiq import h5bse


def test_save(request):
    org_dir = os.getcwd()
    os.chdir(request.fspath.dirname)

    bse = h5bse.h5BSE("test.h5")

    assert bse.save("beta", 0.5)
    assert bse.get("beta") == 0.5

    #Overwrite by another value
    assert bse.save("beta", 0.75)
    assert bse.get("beta") == 0.75

    #Check tuple
    assert bse.save( ("beta") , 0.85)
    assert bse.get("beta") == 0.85


def test_save_data(request):
    org_dir = os.getcwd()
    os.chdir(request.fspath.dirname)

    bse = h5bse.h5BSE("test.h5")

    block_name = [b"up-up", b"up-down", b"down-up", b"down-down"]
    assert bse.save("block_name", block_name)
    assert bse.get("block_name") == block_name

    inner_name = [b"0-0"]
    assert bse.save("inner_name", inner_name)
    assert bse.get("inner_name") == inner_name

    Nin = 4
    Nw = 10
    testdata = numpy.ones((Nin, Nin, Nw, Nw))
    data_loc = {(0, 0): testdata,
                (0, 1): testdata,
                (1, 0): testdata,
                (1, 1): 2.0 * testdata,
                (2, 2): 3.0 * testdata,
                (3, 3): 4.0 * testdata
                }
    X_loc = ("X_loc", 1)
    assert bse.save(X_loc, data_loc)
    X_loc_ref = bse.get(X_loc)
    for _key in list(data_loc.keys()):
        numpy.testing.assert_array_equal(X_loc_ref[_key], data_loc[_key])

    X_loc = ("X_loc", 1, "0.00.00")  # wrong key
    assert bse.save(X_loc, data_loc) == False

    X0q = ("X0_q", 1, "0.00.00")
    testdata = numpy.ones((Nin, Nin, Nw))
    data_C = {(0, 0): testdata,
                (0, 1): testdata,
                (1, 0): testdata,
                (1, 1): 2.0 * testdata,
                (2, 2): 3.0 * testdata,
                (3, 3): 4.0 * testdata
                }
    assert bse.save(X0q, data_C)
    #key_list=bse.get_keylist_data()
    #print(key_list)


def test_save_data_gamma0(request):
    org_dir = os.getcwd()
    os.chdir(request.fspath.dirname)

    bse = h5bse.h5BSE("test.h5")
    Nin = 4
    gamma0 = ("gamma0")
    testdata = numpy.ones((Nin, Nin))
    data_gamma = {(0, 0): testdata,
                (0, 1): testdata,
                (1, 0): testdata,
                (1, 1): 2.0*testdata,
                (2, 2): 3.0*testdata,
                (3, 3): 4.0*testdata
                }
    assert bse.save(gamma0, data_gamma)
    gamma0_ref = bse.get(gamma0)
    for _key in list(data_gamma.keys()):
        numpy.testing.assert_array_equal(gamma0_ref[_key], data_gamma[_key])
