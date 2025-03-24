import numpy
from itertools import product
from collections import Counter
from math import sqrt
from chiq.point_group import PointGroup


# column vectors
cf_d = {
    'E': numpy.array([
        [0, 0, 1, 0, 0],  # 3z^2 - r^2
        [1/sqrt(2), 0, 0, 0, 1/sqrt(2)],  # x^2 - y^2
    ]).T,
    'T2': numpy.array([
        [1/sqrt(2), 0, 0, 0, -1/sqrt(2)],  # xy
        [0, 1/sqrt(2), 0, 1/sqrt(2), 0],  # yz
        [0, 1/sqrt(2), 0, -1/sqrt(2), 0],  # zx
    ]).T
}

cf_f = {
    'E5/2': numpy.array([
        [-sqrt(1/6.), 0, 0, 0, sqrt(5/6.), 0],  # |u>
        [0, sqrt(5/6.), 0, 0, 0, -sqrt(1/6.)],  # |d>
    ]).T,
    'G3/2': numpy.array([
        [sqrt(5/6.), 0, 0, 0, sqrt(1/6.), 0],  # |+u>
        [0, sqrt(1/6.), 0, 0, 0, sqrt(5/6.)],  # |+d>
        [0, 0, 1, 0, 0, 0],  # |-u>
        [0, 0, 0, 1, 0, 0],  # |-d>
    ]).T,
}

cf_refs = {'2': cf_d, '5/2': cf_f}


def _test_cf_states(angular_momentum):
    group = 'O'
    print("\n=========== J=%s group=%s =============" % (repr(angular_momentum), repr(group)))

    pg = PointGroup(angular_momentum=angular_momentum, group=group)
    # pg.G.use_single_reps()

    cf_states = pg.cf_states
    cf_ref = cf_refs[angular_momentum]

    for cf in cf_states:
        print("\n------- cf.name=%s, cf.dim=%d" % (cf.name, cf.dim))
        # print(cf.vecs)
        assert cf.name in cf_ref

        for j in range(cf.dim):
            print("\n#%d" % j)
            vec = cf.vecs[:, j]  # column vector
            print(vec)

            # overlap with reference vectors
            proj = vec.conj() @ cf_ref[cf.name]
            print("projection : <vec|vec_ref> =", proj)
            ol = numpy.linalg.norm(proj)
            print("overlap    : |<vec|vec_ref>|^2 =", ol)
            assert numpy.allclose(ol, 1)

def test_cf_states():
    _test_cf_states(angular_momentum='2')
    _test_cf_states(angular_momentum='5/2')
