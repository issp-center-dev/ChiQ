import numpy as np
from chiq.solver import layout

def _blk(v):
    return np.array([[v]], dtype=complex)

def test_sub_union_of_keys():
    a = {(0, 0): _blk(3), (0, 1): _blk(2)}
    b = {(0, 0): _blk(1), (1, 1): _blk(5)}
    r = layout.sub(a, b)
    assert set(r) == {(0, 0), (0, 1), (1, 1)}
    assert r[(0, 0)] == _blk(2)
    assert r[(0, 1)] == _blk(2)
    assert r[(1, 1)] == _blk(-5)

def test_add_and_scale():
    a = {(0, 0): _blk(3)}
    b = {(0, 0): _blk(1)}
    assert layout.add(a, b)[(0, 0)] == _blk(4)
    assert layout.scale(2.0, a)[(0, 0)] == _blk(6)

def test_identity():
    ident = layout.identity([2, 1])
    assert set(ident) == {(0, 0), (1, 1)}
    assert np.array_equal(ident[(0, 0)], np.eye(2, dtype=complex))
    assert np.array_equal(ident[(1, 1)], np.eye(1, dtype=complex))
