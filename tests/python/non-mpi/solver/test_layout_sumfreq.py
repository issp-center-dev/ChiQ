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
