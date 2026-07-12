import numpy as np
from chiq.solver.kernels import calc_rpa, calc_rrpa

def test_rpa_zero_gamma_gives_chi0():
    nb, nw, n_in, beta = 1, 1, 1, 1.0
    chi0_q = {(0, 0): np.array([[2.0]], complex)}
    gamma0 = {(0, 0): np.array([[0.0]], complex)}
    out = calc_rpa({"chi0_q": chi0_q, "gamma0": gamma0}, beta, nb, nw, n_in)
    assert np.allclose(out["chi_q_rpa"][(0, 0)], [[2.0]])

def test_rrpa_when_chi0loc_equals_chiloc_gives_chi0q():
    # U_eff = chi0_loc^-1 - chi_loc^-1 = 0 -> chi_q_rrpa = chi0_q
    nb, nw, n_in, beta = 1, 1, 1, 1.0
    same = {(0, 0): np.array([[3.0]], complex)}
    chi0_q = {(0, 0): np.array([[2.0]], complex)}
    out = calc_rrpa({"chi0_loc": same, "chi_loc": same, "chi0_q": chi0_q}, beta, nb, nw, n_in)
    assert np.allclose(out["chi_q_rrpa"][(0, 0)], [[2.0]])
