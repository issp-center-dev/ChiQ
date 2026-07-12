"""Block-matrix algebra for the ChiQ solver backends.

Reproduces the A/B/C block conventions and operations of the C++ solver
(src/bse.hpp, src/block_matrix.hpp) on dict-of-ndarray block matrices:
    block matrix := dict[(vi, vj) -> np.ndarray(complex128)]
Missing key == structural zero. See the design spec section 3.4-3.5.
"""

import numpy as np


def parse_matrix_info(info_A, info_B, info_C):
    """Return (nb, nw, n_in) from the three matrix-info vectors.

    Enforces the uniform-dimension invariant (same n_in and nw for every
    block). Raises ValueError on any non-uniform or inconsistent input.
    """
    info_A = list(info_A)
    info_B = list(info_B)
    info_C = list(info_C)
    if not info_C:
        raise ValueError("info_C is empty")
    nb = len(info_C)
    n_in = int(info_C[0])
    if n_in <= 0:
        raise ValueError(f"non-positive inner dimension {n_in}")
    if any(int(x) != n_in for x in info_C):
        raise ValueError("info_C is not uniform")
    if len(info_B) % nb != 0:
        raise ValueError("len(info_B) is not a multiple of nb")
    nw = len(info_B) // nb
    if nw <= 0:
        raise ValueError("nw must be positive")
    if any(int(x) != n_in for x in info_B):
        raise ValueError("info_B is not uniform (!= n_in)")
    if len(info_A) != nb or any(int(x) != n_in * nw for x in info_A):
        raise ValueError("info_A inconsistent with n_in*nw")
    return nb, nw, n_in
