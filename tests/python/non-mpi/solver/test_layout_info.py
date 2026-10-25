import numpy as np
import pytest
from chiq.solver import layout

def test_parse_matrix_info_uniform():
    nb, nw, n_in = 2, 3, 4
    info_A = [n_in * nw] * nb
    info_B = [n_in] * (nb * nw)
    info_C = [n_in] * nb
    assert layout.parse_matrix_info(info_A, info_B, info_C) == (nb, nw, n_in)

def test_parse_matrix_info_rejects_nonuniform_inner():
    with pytest.raises(ValueError):
        layout.parse_matrix_info([8, 8], [4, 4, 5, 4], [4, 4])  # info_B not all equal

def test_parse_matrix_info_rejects_bad_A_dim():
    # nb=2, nw=2, n_in=4 -> info_A must be [8,8]; give [9,9]
    with pytest.raises(ValueError):
        layout.parse_matrix_info([9, 9], [4, 4, 4, 4], [4, 4])
