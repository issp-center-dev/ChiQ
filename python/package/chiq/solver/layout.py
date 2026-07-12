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


def add(a, b):
    """Entrywise a + b over the union of present keys."""
    out = {}
    for k in set(a) | set(b):
        if k in a and k in b:
            out[k] = a[k] + b[k]
        elif k in a:
            out[k] = a[k].copy()
        else:
            out[k] = b[k].copy()
    return out


def sub(a, b):
    """Entrywise a - b over the union of present keys."""
    out = {}
    for k in set(a) | set(b):
        if k in a and k in b:
            out[k] = a[k] - b[k]
        elif k in a:
            out[k] = a[k].copy()
        else:
            out[k] = -b[k]
    return out


def scale(s, bm):
    """Scalar multiply every block by s."""
    return {k: s * v for k, v in bm.items()}


def identity(dims):
    """Block-diagonal identity with the given per-vertex dimensions."""
    return {(i, i): np.eye(int(d), dtype=complex) for i, d in enumerate(dims)}


def connected_components(bm, nvert):
    """Connected components of the block-sparsity graph.

    Edge between i and j iff (i,j) or (j,i) is a present key (matches
    block_matrix.hpp `connected`). Returns a list of components, each a
    sorted list of vertices, components sorted by their first vertex.
    """
    adj = {i: set() for i in range(nvert)}
    for (i, j) in bm:
        adj[i].add(j)
        adj[j].add(i)
    seen = [False] * nvert
    comps = []
    for start in range(nvert):
        if seen[start]:
            continue
        stack = [start]
        seen[start] = True
        comp = []
        while stack:
            v = stack.pop()
            comp.append(v)
            for w in adj[v]:
                if not seen[w]:
                    seen[w] = True
                    stack.append(w)
        comps.append(sorted(comp))
    comps.sort(key=lambda c: c[0])
    return comps


def block_inverse(bm, dims):
    """Per-connected-component dense inverse.

    Assembles each component into a dense super-block, inverts it, and
    scatters back every intra-component vertex pair (dense fill-in,
    numerical zeros retained), matching block_matrix.hpp `inverse`.
    """
    dims = [int(d) for d in dims]
    nvert = len(dims)
    out = {}
    for comp in connected_components(bm, nvert):
        offsets = {}
        size = 0
        for v in comp:
            offsets[v] = size
            size += dims[v]
        super_block = np.zeros((size, size), dtype=complex)
        for a in comp:
            for b in comp:
                if (a, b) in bm:
                    super_block[offsets[a]:offsets[a] + dims[a],
                                offsets[b]:offsets[b] + dims[b]] = bm[(a, b)]
        with np.errstate(divide="ignore", invalid="ignore"):
            try:
                inv = np.linalg.inv(super_block)
            except np.linalg.LinAlgError:
                inv = np.full((size, size), np.nan, dtype=complex)
        for a in comp:
            for b in comp:
                out[(a, b)] = inv[offsets[a]:offsets[a] + dims[a],
                                  offsets[b]:offsets[b] + dims[b]].copy()
    return out
