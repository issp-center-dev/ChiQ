import numpy as np
from chiq.solver import layout

def test_a2b_b2a_roundtrip():
    nb, nw, n_in = 1, 2, 1
    A = np.array([[1, 2], [3, 4]], dtype=complex)  # in*nw+iw indexing, n_in=1
    a = {(0, 0): A}
    b = layout.convert_a2b(a, nb, nw, n_in)
    # keys are (iw1*nb+0, iw2*nb+0) for nonzero extracted entries
    assert (0, 0) in b and (1, 1) in b  # A[0,0]=1, A[1,1]=4
    back = layout.convert_b2a(b, nb, nw, n_in)
    assert np.allclose(back[(0, 0)], A)

def test_a2b_prunes_exact_zero_blocks():
    nb, nw, n_in = 1, 2, 1
    A = np.array([[1, 0], [0, 4]], dtype=complex)
    b = layout.convert_a2b({(0, 0): A}, nb, nw, n_in)
    assert set(b) == {(0, 0), (1, 1)}  # cross-freq (0,1),(1,0) are exactly zero -> pruned

def test_prod_ab_matches_dense():
    nb, nw, n_in = 1, 2, 1
    A = np.array([[1, 2], [3, 4]], dtype=complex)
    # B-type diagonal in iw: keys (iw,iw)
    Bmat = {(0, 0): np.array([[2]], dtype=complex), (1, 1): np.array([[5]], dtype=complex)}
    got = layout.prod_ab(A_as := {(0, 0): A}, Bmat, nb, nw, n_in)
    # dense: convert A to B (full 2x2 in freq), multiply by diag(2,5), convert back
    Bdense = np.diag([2.0, 5.0]).astype(complex)
    assert np.allclose(got[(0, 0)], A @ Bdense)

def test_convert_b2a_prunes_present_but_zero_block():
    # a B key that is present but exactly zero must NOT produce an A block
    # (C++ Convert_B2A prunes on the assembled value via isZero(0), not on key presence)
    B = {(0, 0): np.zeros((1, 1), dtype=complex)}
    out = layout.convert_b2a(B, nb=1, nw=2, n_in=1)
    assert out == {}
