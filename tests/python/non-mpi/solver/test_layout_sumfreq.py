import numpy as np
from chiq.solver import layout

def test_sumfreq_b_sums_over_frequency_and_emits_all_keys():
    nb, nw, n_in = 2, 2, 1
    bm = {}
    # only diagonal (i==j) same-iw keys present
    bm[(0 * nb + 0, 0 * nb + 0)] = np.array([[1.0]], dtype=complex)
    bm[(1 * nb + 0, 1 * nb + 0)] = np.array([[3.0]], dtype=complex)  # iw=1, b=0
    out = layout.sumfreq_b(bm, nb, nw, n_in)
    assert set(out) == {(i, j) for i in range(nb) for j in range(nb)}
    assert np.allclose(out[(0, 0)], [[4.0]])  # 1 + 3
    assert np.allclose(out[(0, 1)], [[0.0]])

def test_sumfreq_a_sums_over_both_frequency_indices():
    nb, nw, n_in = 1, 2, 1
    # A block is (n_in*nw)=(2) x 2; inner flattened as in*nw+iw -> here in=0 so index==iw
    A = np.array([[1, 2], [3, 4]], dtype=complex)
    bm = {(0, 0): A}
    out = layout.sumfreq_a(bm, nb, nw, n_in)
    assert set(out) == {(0, 0)}
    assert np.allclose(out[(0, 0)], [[1 + 2 + 3 + 4]])

def test_sumfreq_a_multiorbital_reshape_order():
    # n_in=2, nw=2 -> A block is 4x4; distinguishes a*nw+iw from the transposed iw*n_in+a.
    nb, nw, n_in = 1, 2, 2
    A = np.arange(16, dtype=complex).reshape(4, 4)
    bm = {(0, 0): A}
    out = layout.sumfreq_a(bm, nb, nw, n_in)
    # explicit reference: expected[a,c] = sum over iw1,iw2 of A[a*nw+iw1, c*nw+iw2]
    expected = np.zeros((n_in, n_in), dtype=complex)
    for a in range(n_in):
        for c in range(n_in):
            for iw1 in range(nw):
                for iw2 in range(nw):
                    expected[a, c] += A[a * nw + iw1, c * nw + iw2]
    assert np.allclose(out[(0, 0)], expected)

def test_sumfreq_b_offdiagonal_multiorbital():
    # n_in=2, nb=2, nw=2: off-diagonal (i!=j) blocks summed over iw.
    nb, nw, n_in = 2, 2, 2
    b0 = np.arange(4, dtype=complex).reshape(2, 2)
    b1 = 10 + np.arange(4, dtype=complex).reshape(2, 2)
    # key (iw*nb + i, iw*nb + j) with i=0, j=1
    bm = {(0 * nb + 0, 0 * nb + 1): b0, (1 * nb + 0, 1 * nb + 1): b1}
    out = layout.sumfreq_b(bm, nb, nw, n_in)
    assert set(out) == {(i, j) for i in range(nb) for j in range(nb)}
    assert np.allclose(out[(0, 1)], b0 + b1)
    assert np.allclose(out[(0, 0)], np.zeros((n_in, n_in)))
