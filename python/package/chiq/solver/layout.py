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


def matmul(a, b, dims):
    """Block matrix product a @ b.

    Emits key (i,j) iff there exists k with (i,k) in a and (k,j) in b
    (structural rule, matching block_matrix.hpp operator*); the numerical
    value never decides key presence.
    """
    nvert = len(dims)
    out = {}
    for i in range(nvert):
        for j in range(nvert):
            acc = None
            for k in range(nvert):
                if (i, k) in a and (k, j) in b:
                    p = a[(i, k)] @ b[(k, j)]
                    acc = p if acc is None else acc + p
            if acc is not None:
                out[(i, j)] = acc
    return out


def sumfreq_a(bm, nb, nw, n_in):
    """A-type -> C-type by summing over both frequency indices.

    out[(i,j)][a,c] = sum_{iw1,iw2} A[(i,j)][a*nw+iw1, c*nw+iw2].
    Only keys present in bm are emitted (matches Sumfreq_A which skips
    absent A blocks).
    """
    out = {}
    for (i, j), mat in bm.items():
        r = mat.reshape(n_in, nw, n_in, nw)
        out[(i, j)] = r.sum(axis=(1, 3))
    return out


def sumfreq_b(bm, nb, nw, n_in):
    """B-type -> C-type by summing over the (diagonal) frequency index.

    out[(i,j)] = sum_iw B[(iw*nb+i, iw*nb+j)] over present keys. Every
    nb*nb key is emitted (matches Sumfreq_B, which assigns all blocks).
    """
    out = {}
    for i in range(nb):
        for j in range(nb):
            m = np.zeros((n_in, n_in), dtype=complex)
            for iw in range(nw):
                key = (iw * nb + i, iw * nb + j)
                if key in bm:
                    m = m + bm[key]
            out[(i, j)] = m
    return out


def convert_a2b(bm, nb, nw, n_in):
    """A-type -> B-type. Emits (iw1*nb+b1, iw2*nb+b2) iff the extracted
    block is not exactly zero (matches Convert_A2B isZero(0) pruning)."""
    out = {}
    for (b1, b2), mat in bm.items():
        r = mat.reshape(n_in, nw, n_in, nw)
        for iw1 in range(nw):
            for iw2 in range(nw):
                block = r[:, iw1, :, iw2]
                if np.any(block != 0):
                    out[(iw1 * nb + b1, iw2 * nb + b2)] = block.copy()
    return out


def convert_b2a(bm, nb, nw, n_in):
    """B-type -> A-type. Emits (b1,b2) iff its assembled A block is not
    exactly zero (matches Convert_B2A isZero(0) pruning)."""
    out = {}
    for b1 in range(nb):
        for b2 in range(nb):
            r = np.zeros((n_in, nw, n_in, nw), dtype=complex)
            for iw1 in range(nw):
                for iw2 in range(nw):
                    key = (iw1 * nb + b1, iw2 * nb + b2)
                    if key in bm:
                        r[:, iw1, :, iw2] = bm[key]
            if np.any(r != 0):
                out[(b1, b2)] = r.reshape(n_in * nw, n_in * nw)
    return out


def prod_ab(a_mat, b_mat, nb, nw, n_in):
    """A-type * B-type -> A-type, via convert to B, multiply, convert back
    (reproduces bse.hpp Prod_A_B)."""
    a_in_b = convert_a2b(a_mat, nb, nw, n_in)
    ab_in_b = matmul(a_in_b, b_mat, [n_in] * (nb * nw))
    return convert_b2a(ab_in_b, nb, nw, n_in)
