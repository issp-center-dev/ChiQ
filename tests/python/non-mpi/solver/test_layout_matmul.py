import numpy as np
from chiq.solver import layout

def test_matmul_structural_keys_and_values():
    a = {(0, 0): np.array([[1, 2], [3, 4]], dtype=complex),
         (0, 1): np.array([[1, 0], [0, 1]], dtype=complex)}
    b = {(0, 0): np.array([[1, 0], [0, 1]], dtype=complex),
         (1, 1): np.array([[2, 0], [0, 2]], dtype=complex)}
    r = layout.matmul(a, b, dims=[2, 2])
    # (0,0): a00@b00 ; (0,1): a01@b11 ; nothing writes (1,*)
    assert set(r) == {(0, 0), (0, 1)}
    assert np.allclose(r[(0, 0)], a[(0, 0)])
    assert np.allclose(r[(0, 1)], 2 * np.eye(2))

def test_matmul_noncommuting():
    x = np.array([[0, 1], [0, 0]], dtype=complex)
    y = np.array([[0, 0], [1, 0]], dtype=complex)
    a = {(0, 0): x}
    b = {(0, 0): y}
    assert np.allclose(layout.matmul(a, b, [2])[(0, 0)], x @ y)
    assert not np.allclose(x @ y, y @ x)
