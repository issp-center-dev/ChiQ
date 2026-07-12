import numpy as np
from chiq.solver import get_solver

def test_chi0_formula_one_block():
    # nb=1, nw=1, n_in=1, beta given; chi0_q = chi0_loc + (1/beta)*(X0_q - X0_loc)
    beta = 4.0
    iA, iB, iC = [1], [1], [1]
    s = get_solver("numpy", beta, iA, iB, iC)
    s.set({(0, 0): np.array([[7.0]], complex)}, "chi0_loc")
    s.set({(0, 0): np.array([[2.0]], complex)}, "X0_loc")
    s.set({(0, 0): np.array([[10.0]], complex)}, "X0_q")
    s.calc("chi0")
    out = s.get("chi0_q")
    assert np.allclose(out[(0, 0)], [[7.0 + (10.0 - 2.0) / beta]])
