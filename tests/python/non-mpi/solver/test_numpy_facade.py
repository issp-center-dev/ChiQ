import numpy as np
import pytest
from chiq.solver import get_solver

def _info(nb, nw, n_in):
    return ([n_in * nw] * nb, [n_in] * (nb * nw), [n_in] * nb)

def test_set_unknown_name_raises():
    s = get_solver("numpy", 10.0, *_info(1, 2, 1))
    with pytest.raises(RuntimeError, match="Invalid type"):
        s.set({(0, 0): np.zeros((1, 1), complex)}, "nope")

def test_get_before_calc_returns_empty_dict():
    s = get_solver("numpy", 10.0, *_info(1, 2, 1))
    assert s.get("chi0_q") == {}

def test_missing_input_reports_capital_X0_Loc():
    s = get_solver("numpy", 10.0, *_info(1, 2, 1))
    s.set({(0, 0): np.zeros((1, 1), complex)}, "chi0_loc")
    s.set({(0, 0): np.zeros((1, 1), complex), (1, 1): np.zeros((1, 1), complex)}, "X0_q")
    with pytest.raises(RuntimeError, match="'X0_Loc' must be set"):
        s.calc("chi0")

def test_missing_phi_sum_reports_capital_Phi_Sum():
    # SCL requires X0_q, X0_loc, Phi, Phi_sum; omit Phi_sum -> error names 'Phi_Sum' (capital S), matching C++.
    nb, nw, n_in = 1, 2, 1
    s = get_solver("numpy", 10.0, *_info(nb, nw, n_in))
    X0 = {(0, 0): np.zeros((n_in, n_in), complex), (1, 1): np.zeros((n_in, n_in), complex)}
    s.set(X0, "X0_q")
    s.set(X0, "X0_loc")
    s.set(X0, "Phi")
    with pytest.raises(RuntimeError, match="'Phi_Sum' must be set"):
        s.calc("scl")

def test_get_unknown_name_raises():
    s = get_solver("numpy", 10.0, *_info(1, 2, 1))
    with pytest.raises(RuntimeError, match="Invalid type"):
        s.get("totally_bogus_name")

def test_get_valid_uncomputed_name_returns_empty_dict():
    s = get_solver("numpy", 10.0, *_info(1, 2, 1))
    assert s.get("chi_q") == {}   # valid output name, not computed yet -> {}

def test_set_copies_input_no_aliasing():
    # mutating the caller's array after set() must not change the stored input / results
    nb, nw, n_in, beta = 1, 1, 1, 4.0
    s = get_solver("numpy", beta, [1], [1], [1])
    chi0_loc = {(0, 0): np.array([[7.0]], complex)}
    X0_loc = {(0, 0): np.array([[2.0]], complex)}
    X0_q = {(0, 0): np.array([[10.0]], complex)}
    s.set(chi0_loc, "chi0_loc"); s.set(X0_loc, "X0_loc"); s.set(X0_q, "X0_q")
    X0_q[(0, 0)][0, 0] = 999.0   # mutate caller's array AFTER set
    s.calc("chi0")
    out = s.get("chi0_q")
    # result uses the value at set() time (10.0), not the mutated 999.0
    assert np.allclose(out[(0, 0)], [[7.0 + (10.0 - 2.0) / beta]])

def test_get_returns_owned_copy_no_aliasing():
    nb, nw, n_in, beta = 1, 1, 1, 4.0
    s = get_solver("numpy", beta, [1], [1], [1])
    s.set({(0, 0): np.array([[7.0]], complex)}, "chi0_loc")
    s.set({(0, 0): np.array([[2.0]], complex)}, "X0_loc")
    s.set({(0, 0): np.array([[10.0]], complex)}, "X0_q")
    s.calc("chi0")
    out1 = s.get("chi0_q")
    out1[(0, 0)][0, 0] = -12345.0   # mutate the returned array
    out2 = s.get("chi0_q")
    assert np.allclose(out2[(0, 0)], [[7.0 + (10.0 - 2.0) / beta]])  # unchanged
