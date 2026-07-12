import numpy as np
from chiq.solver import get_solver
from chiq.solver.kernels import calc_bse

def _info(nb, nw, n_in):
    return ([n_in * nw] * nb, [n_in] * (nb * nw), [n_in] * nb)

def test_iq_and_iq_scl_share_one_slot():
    # After a bse calc, get("I_q_scl") returns the same value as get("I_q").
    nb, nw, n_in = 1, 2, 1
    s = get_solver("numpy", 5.0, *_info(nb, nw, n_in))
    X0 = {(0, 0): np.array([[1.0]], complex), (1, 1): np.array([[2.0]], complex)}
    s.set(X0, "X0_q")
    s.set(X0, "X0_loc")
    s.set({(0, 0): np.eye(2, dtype=complex)}, "X_loc")
    s.set({(0, 0): np.array([[3.0]], complex)}, "chi_loc")
    s.calc("bse")
    assert s.get("I_q_scl") == s.get("I_q")

def test_bse_reduces_to_chi_loc_when_X0q_equals_X0loc():
    # If X0_q == X0_loc then P_q = 0, Z_q = 0, so chi_q == chi_loc.
    nb, nw, n_in, beta = 1, 2, 1, 5.0
    X0 = {(0, 0): np.array([[1.0]], complex), (1, 1): np.array([[2.0]], complex)}
    X_loc = {(0, 0): np.array([[1.0, 0.0], [0.0, 1.0]], complex)}  # A-type 2x2
    chi_loc = {(0, 0): np.array([[3.0]], complex)}
    state = {"X0_q": X0, "X0_loc": X0, "X_loc": X_loc, "chi_loc": chi_loc}
    out = calc_bse(state, beta, nb, nw, n_in)
    assert np.allclose(out["chi_q"][(0, 0)], [[3.0]])
    # I_q = chi_loc^-1 - chi_q^-1 = 0 here
    assert np.allclose(out["I_q"][(0, 0)], [[0.0]])
