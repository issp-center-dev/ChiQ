import numpy as np
from chiq.solver import layout

def test_connected_components_two_islands():
    bm = {(0, 0): np.eye(1) * 1, (2, 2): np.eye(1) * 1, (0, 2): np.eye(1) * 0.0 + 1}
    # vertices 0 and 2 are connected via (0,2); vertex 1 isolated
    comps = layout.connected_components(bm, nvert=3)
    assert comps == [[0, 2], [1]]

def test_block_inverse_matches_dense_and_fills_in():
    # single 2-vertex component; inverse is generally full within the component
    a = np.array([[2, 1], [0, 3]], dtype=complex)
    b = np.array([[1, 0], [1, 1]], dtype=complex)
    c = np.array([[4, 0], [0, 2]], dtype=complex)
    d = np.array([[1, 2], [0, 1]], dtype=complex)
    bm = {(0, 0): a, (0, 1): b, (1, 0): c, (1, 1): d}
    inv = layout.block_inverse(bm, dims=[2, 2])
    # reconstruct dense 4x4 and compare
    big = np.block([[a, b], [c, d]])
    big_inv = np.linalg.inv(big)
    got = np.block([[inv[(0, 0)], inv[(0, 1)]], [inv[(1, 0)], inv[(1, 1)]]])
    assert np.allclose(got, big_inv)
    # fill-in: all 4 intra-component pairs present
    assert set(inv) == {(0, 0), (0, 1), (1, 0), (1, 1)}

def test_block_inverse_keeps_components_separate():
    bm = {(0, 0): np.array([[2]], dtype=complex), (1, 1): np.array([[4]], dtype=complex)}
    inv = layout.block_inverse(bm, dims=[1, 1])
    assert set(inv) == {(0, 0), (1, 1)}
    assert np.allclose(inv[(0, 0)], [[0.5]])
    assert np.allclose(inv[(1, 1)], [[0.25]])

def test_block_inverse_singular_gives_nan():
    bm = {(0, 0): np.zeros((2, 2), dtype=complex)}
    inv = layout.block_inverse(bm, dims=[2])
    assert set(inv) == {(0, 0)}
    assert np.isnan(inv[(0, 0)]).all()

def test_block_inverse_fills_absent_intra_component_key():
    a = np.array([[2.0]], dtype=complex)
    c = np.array([[1.0]], dtype=complex)   # (1,0) present -> vertices 0,1 connected
    d = np.array([[3.0]], dtype=complex)
    bm = {(0, 0): a, (1, 0): c, (1, 1): d}   # (0,1) ABSENT in input
    inv = layout.block_inverse(bm, dims=[1, 1])
    assert (0, 1) in inv   # fill-in added the absent key
    # compare to dense inverse treating missing (0,1) as zero
    big = np.block([[a, np.zeros((1, 1), complex)], [c, d]])
    big_inv = np.linalg.inv(big)
    got = np.block([[inv[(0, 0)], inv[(0, 1)]], [inv[(1, 0)], inv[(1, 1)]]])
    assert np.allclose(got, big_inv)
