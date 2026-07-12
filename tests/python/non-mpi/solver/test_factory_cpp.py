import numpy as np
import pytest
from chiq.solver import get_solver

def _info(nb, nw, n_in):
    return ([n_in * nw] * nb, [n_in] * (nb * nw), [n_in] * nb)

def test_factory_returns_cpp_by_default_and_runs_chi0():
    nb, nw, n_in = 1, 2, 1
    iA, iB, iC = _info(nb, nw, n_in)
    solver = get_solver("cpp", beta=10.0, info_A=iA, info_B=iB, info_C=iC)
    # minimal chi0 inputs
    solver.set({(0, 0): np.zeros((n_in, n_in), dtype=complex)}, "chi0_loc")
    solver.set({(0, 0): np.zeros((n_in, n_in), dtype=complex), (1, 1): np.zeros((n_in, n_in), dtype=complex)}, "X0_loc")
    solver.set({(0, 0): np.zeros((n_in, n_in), dtype=complex), (1, 1): np.zeros((n_in, n_in), dtype=complex)}, "X0_q")
    solver.calc("chi0")
    out = solver.get("chi0_q")
    assert isinstance(out, dict)

def test_factory_rejects_unknown_backend():
    iA, iB, iC = _info(1, 2, 1)
    with pytest.raises(ValueError):
        get_solver("gpu999", 10.0, iA, iB, iC)
