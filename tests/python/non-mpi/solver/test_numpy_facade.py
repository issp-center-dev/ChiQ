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

def test_get_unknown_name_raises():
    s = get_solver("numpy", 10.0, *_info(1, 2, 1))
    with pytest.raises(RuntimeError, match="Invalid type"):
        s.get("totally_bogus_name")

def test_get_valid_uncomputed_name_returns_empty_dict():
    s = get_solver("numpy", 10.0, *_info(1, 2, 1))
    assert s.get("chi_q") == {}   # valid output name, not computed yet -> {}
