import numpy as np
from chiq.solver.kernels import calc_scl

def test_scl_when_X0q_equals_X0loc_gives_phisum_squared():
    # Lambda_q = X0_loc^-1 - X0_q^-1 = 0 -> lambda_q = 0
    # chi_q_scl = Phi_sum * I^-1 * Phi_sum = Phi_sum @ Phi_sum
    nb, nw, n_in, beta = 1, 2, 1, 1.0
    X0 = {(0, 0): np.array([[1.0]], complex), (1, 1): np.array([[2.0]], complex)}
    Phi = {(0, 0): np.array([[1.0]], complex), (1, 1): np.array([[1.0]], complex)}  # B-type
    Phi_sum = {(0, 0): np.array([[3.0]], complex)}  # C-type
    out = calc_scl({"X0_q": X0, "X0_loc": X0, "Phi": Phi, "Phi_sum": Phi_sum}, beta, nb, nw, n_in)
    assert np.allclose(out["chi_q_scl"][(0, 0)], [[9.0]])   # 3*3
    assert np.allclose(out["I_q_scl"][(0, 0)], [[0.0]])     # Phi_sum^-1 * 0 * Phi_sum^-1
