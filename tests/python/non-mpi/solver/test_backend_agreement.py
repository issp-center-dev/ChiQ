import numpy as np
import pytest
from chiq.solver import get_solver

RNG = np.random.default_rng(0)

def _info(nb, nw, n_in):
    return ([n_in * nw] * nb, [n_in] * (nb * nw), [n_in] * nb)

def _rand_c(shape):
    return (RNG.standard_normal(shape) + 1j * RNG.standard_normal(shape))

def _well_conditioned(n):
    # diagonally dominant -> invertible, bounded condition number
    m = _rand_c((n, n))
    return m + n * np.eye(n)

def _assert_agree(a, b):
    for k in set(a) | set(b):
        if k in a and k in b:
            va, vb = a[k], b[k]
        elif k in a:
            va, vb = a[k], np.zeros_like(a[k])
        else:
            va, vb = np.zeros_like(b[k]), b[k]
        # NaN / signed-inf position masks, then elementwise mixed tolerance.
        # np.isposinf/np.isneginf reject complex dtype outright (ambiguous
        # sign), so check real/imag parts separately -- this is a test-harness
        # fix, not a loosening of the tolerance/metric.
        va_r, va_i = np.real(va), np.imag(va)
        vb_r, vb_i = np.real(vb), np.imag(vb)
        assert np.array_equal(np.isnan(va_r), np.isnan(vb_r))
        assert np.array_equal(np.isnan(va_i), np.isnan(vb_i))
        assert np.array_equal(np.isposinf(va_r), np.isposinf(vb_r))
        assert np.array_equal(np.isposinf(va_i), np.isposinf(vb_i))
        assert np.array_equal(np.isneginf(va_r), np.isneginf(vb_r))
        assert np.array_equal(np.isneginf(va_i), np.isneginf(vb_i))
        np.testing.assert_allclose(va, vb, rtol=1e-8, atol=1e-10, equal_nan=True)

def _c_blocks(nb, n_in):
    return {(i, j): _rand_c((n_in, n_in)) for i in range(nb) for j in range(nb)}

def _c_blocks_wc(nb, n_in):
    # C-type with well-conditioned diagonal blocks (invertible)
    d = {(i, i): _well_conditioned(n_in) for i in range(nb)}
    return d

def _b_blocks_diag(nb, nw, n_in):
    return {(iw * nb + b, iw * nb + b2): _rand_c((n_in, n_in))
            for iw in range(nw) for b in range(nb) for b2 in range(nb)}

def _b_blocks_diag_wc(nb, nw, n_in):
    # B-type well-conditioned in each same-iw super-block (diagonal blocks dominant)
    out = {}
    for iw in range(nw):
        for b in range(nb):
            for b2 in range(nb):
                m = _rand_c((n_in, n_in))
                if b == b2:
                    m = m + (nb * n_in) * np.eye(n_in)
                out[(iw * nb + b, iw * nb + b2)] = m
    return out

def _a_blocks(nb, nw, n_in):
    return {(i, j): _rand_c((n_in * nw, n_in * nw)) for i in range(nb) for j in range(nb)}

@pytest.mark.parametrize("nb,nw,n_in", [(1, 2, 1), (2, 2, 2), (2, 3, 2)])
def test_chi0_agreement(nb, nw, n_in):
    iA, iB, iC = _info(nb, nw, n_in)
    inputs = {
        "chi0_loc": {(i, j): _rand_c((n_in, n_in)) for i in range(nb) for j in range(nb)},
        "X0_loc": {(iw * nb + b, iw * nb + b2): _rand_c((n_in, n_in))
                   for iw in range(nw) for b in range(nb) for b2 in range(nb)},
    }
    inputs["X0_q"] = {k: _rand_c((n_in, n_in)) for k in inputs["X0_loc"]}
    outs = {}
    for backend in ("cpp", "numpy"):
        s = get_solver(backend, 12.0, iA, iB, iC)
        for name, bm in inputs.items():
            s.set(bm, name)
        s.calc("chi0")
        outs[backend] = s.get("chi0_q")
    _assert_agree(outs["cpp"], outs["numpy"])

@pytest.mark.parametrize("nb,nw,n_in", [(1, 2, 1), (2, 2, 2)])
def test_rpa_agreement(nb, nw, n_in):
    iA, iB, iC = _info(nb, nw, n_in)
    chi0_q = {(i, j): _rand_c((n_in, n_in)) for i in range(nb) for j in range(nb)}
    gamma0 = {(i, i): 0.01 * _rand_c((n_in, n_in)) for i in range(nb)}
    outs = {}
    for backend in ("cpp", "numpy"):
        s = get_solver(backend, 12.0, iA, iB, iC)
        s.set(chi0_q, "chi0_q")
        s.set(gamma0, "gamma0")
        s.calc("rpa")
        outs[backend] = s.get("chi_q_rpa")
    _assert_agree(outs["cpp"], outs["numpy"])

@pytest.mark.parametrize("nb,nw,n_in", [(1, 2, 1), (2, 2, 2)])
def test_rrpa_agreement(nb, nw, n_in):
    iA, iB, iC = _info(nb, nw, n_in)
    inputs = {
        "chi0_loc": _c_blocks_wc(nb, n_in),
        "chi_loc": _c_blocks_wc(nb, n_in),
        "chi0_q": _c_blocks(nb, n_in),
    }
    outs = {}
    for backend in ("cpp", "numpy"):
        s = get_solver(backend, 12.0, iA, iB, iC)
        for name, bm in inputs.items():
            s.set(bm, name)
        s.calc("rrpa")
        outs[backend] = s.get("chi_q_rrpa")
    _assert_agree(outs["cpp"], outs["numpy"])

@pytest.mark.parametrize("nb,nw,n_in", [(1, 2, 1), (2, 2, 2)])
def test_bse_agreement(nb, nw, n_in):
    iA, iB, iC = _info(nb, nw, n_in)
    X0_loc = _b_blocks_diag_wc(nb, nw, n_in)
    inputs = {
        "X0_loc": X0_loc,
        "X0_q": _b_blocks_diag_wc(nb, nw, n_in),
        "X_loc": _a_blocks(nb, nw, n_in),
        "chi_loc": _c_blocks_wc(nb, n_in),
    }
    outs_chi, outs_iq = {}, {}
    for backend in ("cpp", "numpy"):
        s = get_solver(backend, 12.0, iA, iB, iC)
        for name, bm in inputs.items():
            s.set(bm, name)
        s.calc("bse")
        outs_chi[backend] = s.get("chi_q")
        outs_iq[backend] = s.get("I_q")
    _assert_agree(outs_chi["cpp"], outs_chi["numpy"])
    _assert_agree(outs_iq["cpp"], outs_iq["numpy"])

@pytest.mark.parametrize("nb,nw,n_in", [(1, 2, 1), (2, 2, 2)])
def test_scl_agreement(nb, nw, n_in):
    iA, iB, iC = _info(nb, nw, n_in)
    inputs = {
        "X0_loc": _b_blocks_diag_wc(nb, nw, n_in),
        "X0_q": _b_blocks_diag_wc(nb, nw, n_in),
        "Phi": _b_blocks_diag(nb, nw, n_in),
        "Phi_sum": _c_blocks_wc(nb, n_in),
    }
    outs_chi, outs_iq = {}, {}
    for backend in ("cpp", "numpy"):
        s = get_solver(backend, 12.0, iA, iB, iC)
        for name, bm in inputs.items():
            s.set(bm, name)
        s.calc("scl")
        outs_chi[backend] = s.get("chi_q_scl")
        outs_iq[backend] = s.get("I_q_scl")
    _assert_agree(outs_chi["cpp"], outs_chi["numpy"])
    _assert_agree(outs_iq["cpp"], outs_iq["numpy"])
