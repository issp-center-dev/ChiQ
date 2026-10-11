"""Reconstruct X0(q) from the k-resolved Green's function G(k).

This module does not depend on DCore. It defines the HDF5 layout of the
``x0_source`` group written by ``dcore_chiq`` with ``X0q_mode = "gk"`` and the
reconstruction used by ``chiq_main``.

Convention (matches the default FFT writer ``SumkDFTChi.calc_chi0_fft``)::

    X0[(b_l, b_r)][(i1, i2), (i3, i4), t]
        = -(1/n_k) sum_k G[b1][c3, c1, i3, i1](t, k) * G[b2][c2, c4, i2, i4](t + wb, k - q)

with ``block2_pairs[b_l] = (c1, b1, c2, b2)`` and ``block2_pairs[b_r] = (c3, b3, c4, b4)``;
only pairs with ``b1 == b3`` and ``b2 == b4`` are stored.
"""

import h5py
import numpy as np

FORMAT = "gk"
_KQ_CACHE_SIZE = 16  # k-q permutations kept per process (n_k integers each)
_CHUNK_BYTES = 16 * 2**20  # target size of one G chunk (HDF5 requires < 4 GiB)
VERSION = 1
GROUP = "x0_source"


def _g_chunks(n_iwG, n_k, chunk_bytes=_CHUNK_BYTES):
    """Chunk shape for G/<b>: one (c1, c2, i1, i2) slice, split along k and then iw
    so that a chunk stays below ``chunk_bytes`` for any mesh size."""
    n_elem = max(1, chunk_bytes // np.dtype(np.complex128).itemsize)
    ck = min(n_k, n_elem)
    cw = min(n_iwG, max(1, n_elem // ck))
    return (1, 1, 1, 1, cw, ck)


class X0SourceError(Exception):
    """Invalid or inconsistent X0 input representation."""


def bse_root(groupname=""):
    return "bse/" + groupname if groupname else "bse"


def parse_q_label(label, div):
    """'qx.qy.qz' -> (qx, qy, qz), checked against the k mesh."""
    parts = str(label).split(".")
    if len(parts) != 3:
        raise X0SourceError(f"q label {label!r} is not of the form 'qx.qy.qz'")
    try:
        q = tuple(int(p) for p in parts)
    except ValueError:
        raise X0SourceError(f"q label {label!r} is not of the form 'qx.qy.qz'")
    for ql, nl in zip(q, div):
        if not 0 <= ql < nl:
            raise X0SourceError(f"q label {label!r} is not on the k mesh {tuple(div)}")
    return q


def k_minus_q_index(div, q):
    """Return idx with idx[k] = flat index of (k - q) mod div."""
    kx, ky, kz = np.indices(tuple(div)).reshape(3, -1)
    kx = (kx - q[0]) % div[0]
    ky = (ky - q[1]) % div[1]
    kz = (kz - q[2]) % div[2]
    return kz + div[2] * (ky + div[1] * kx)


def stored_block_pairs(block2_pairs):
    """List of (b_l, b_r) kept in X0: b1 == b3 and b2 == b4."""
    pairs = []
    for b_l, (_, b1, _, b2) in enumerate(block2_pairs):
        for b_r, (_, b3, _, b4) in enumerate(block2_pairs):
            if b1 == b3 and b2 == b4:
                pairs.append((b_l, b_r))
    return pairs


def x0_from_g(G, block2_pairs, inner2_pairs, wb, num_wf, iw_offset, kq_index):
    """Compute X0 for all stored block pairs.

    G: list (by single-particle block) of arrays [c, c', i, i', iw, k]
    kq_index: idx[k] = index of k - q (length n_k)
    Returns {(b_l, b_r): complex array [n_in2, n_in2, 2*num_wf]}.
    """
    n_wf2 = 2 * num_wf
    lo = iw_offset - num_wf
    n_k = len(kq_index)
    a = np.asarray(inner2_pairs)[:, 0]
    b = np.asarray(inner2_pairs)[:, 1]
    result = {}
    for b_l, b_r in stored_block_pairs(block2_pairs):
        c1, b1, c2, b2 = block2_pairs[b_l]
        c3, _, c4, _ = block2_pairs[b_r]
        # gl[i3, i1, t, k], gr[i2, i4, t, k]
        gl = G[b1][c3, c1, :, :, lo:lo + n_wf2, :]
        gr = G[b2][c2, c4, :, :, lo + wb:lo + wb + n_wf2, :][..., kq_index]
        n_in = gl.shape[0]
        # A[t, (i1, i3), k] @ B[t, k, (i2, i4)] -> M[t, (i1, i3), (i2, i4)]
        A = gl.transpose(2, 1, 0, 3).reshape(n_wf2, n_in * n_in, n_k)
        B = gr.transpose(2, 3, 0, 1).reshape(n_wf2, n_k, n_in * n_in)
        M = np.matmul(A, B).reshape(n_wf2, n_in, n_in, n_in, n_in)
        M *= -1.0 / n_k
        # M[t, i1, i3, i2, i4] -> X[(i1, i2), (i3, i4), t]
        x = M[:, a[:, None], a[None, :], b[:, None], b[None, :]]
        result[(b_l, b_r)] = np.ascontiguousarray(x.transpose(1, 2, 0))
    return result


def x0_loc_from_gbar(Gbar, block2_pairs, inner2_pairs, wb, num_wf, iw_offset):
    """X0 averaged over all q of the mesh, from k-averaged G [c, c', i, i', iw]."""
    G1 = [g[..., np.newaxis] for g in Gbar]
    return x0_from_g(G1, block2_pairs, inner2_pairs, wb, num_wf, iw_offset, np.zeros(1, dtype=int))


def _read_name_list(h5f, path):
    grp = h5f[path]
    return [grp[str(i)][()] for i in range(len(grp))]


def detect_representation(h5f, groupname=""):
    """Return 'saved', 'gk', or None (no X0 input at all)."""
    root = bse_root(groupname)
    has_saved = (root + "/input/X0_q") in h5f
    has_gk = (root + "/" + GROUP) in h5f
    if has_saved and has_gk:
        raise X0SourceError(
            "the file contains both 'X0_q' and 'x0_source'. "
            "Rerun dcore_chiq to write a single X0 representation."
        )
    if has_gk:
        attrs = h5f[root + "/" + GROUP].attrs
        if attrs.get("format") != FORMAT or int(attrs.get("version", -1)) != VERSION:
            raise X0SourceError(
                f"unsupported x0_source format {attrs.get('format')!r} version {attrs.get('version')!r}"
            )
        if int(attrs.get("complete", 0)) != 1:
            raise X0SourceError("x0_source is incomplete (dcore_chiq was interrupted). Rerun dcore_chiq.")
        return "gk"
    if has_saved:
        return "saved"
    return None


def remove_x0_source(h5f, groupname=""):
    """Delete an x0_source group (used when X0(q) is stored instead)."""
    path = bse_root(groupname) + "/" + GROUP
    if path in h5f:
        del h5f[path]


def _int_array(grp, name, ncol):
    arr = grp[name][()]
    if arr.ndim != 2 or arr.shape[1] != ncol or not np.issubdtype(arr.dtype, np.integer):
        raise X0SourceError(f"x0_source/{name} must be an integer array of shape (n, {ncol})")
    return arr.astype(int)


class GkX0Source:
    """X0(q) provider backed by the ``x0_source`` group of a gk-mode file."""

    def __init__(self, filename, groupname=""):
        root = bse_root(groupname)
        with h5py.File(filename, "r") as h5f:
            if detect_representation(h5f, groupname) != "gk":
                raise X0SourceError(f"{filename!r} has no x0_source group")
            grp = h5f[root + "/" + GROUP]
            attrs = grp.attrs
            self.num_wf = int(attrs["num_wf"])
            self.list_wb = [int(w) for w in np.atleast_1d(attrs["list_wb"])]
            self.n_iwG = int(attrs["n_iwG"])
            self.iw_offset = int(attrs["iw_offset"])
            self.div = tuple(int(d) for d in attrs["div"])
            beta = float(attrs["beta"])

            n_b1 = len(grp["block1_names"])
            self.block2_pairs = [tuple(p) for p in _int_array(grp, "block2_pairs", 4)]
            self.inner2_pairs = _int_array(grp, "inner2_pairs", 2)
            self.G = []
            for b in range(n_b1):
                if f"G/{b}" not in grp:
                    raise X0SourceError(f"x0_source/G/{b} is missing")
                self.G.append(grp[f"G/{b}"][()])

            n_block = len(h5f[root + "/info/block_name"])
            n_inner = len(h5f[root + "/info/inner_name"])
            info_beta = float(h5f[root + "/info/beta"][()])
            x0_loc_wb = []
            if root + "/input/X0_loc" in h5f:
                x0_loc_wb = list(h5f[root + "/input/X0_loc"].keys())

        self._validate(n_block, n_inner, beta, info_beta, x0_loc_wb)
        self._kq_cache = {}

    def _validate(self, n_block, n_inner, beta, info_beta, x0_loc_wb):
        n_k = int(np.prod(self.div))
        g0 = self.G[0]
        if g0.ndim != 6:
            raise X0SourceError("x0_source/G arrays must have 6 dimensions")
        n_c, n_in = g0.shape[0], g0.shape[2]
        for b, g in enumerate(self.G):
            if g.shape != (n_c, n_c, n_in, n_in, self.n_iwG, n_k) or not np.iscomplexobj(g):
                raise X0SourceError(
                    f"x0_source/G/{b} has shape {g.shape}, expected complex "
                    f"{(n_c, n_c, n_in, n_in, self.n_iwG, n_k)}"
                )
        if len(self.block2_pairs) != n_block:
            raise X0SourceError("x0_source/block2_pairs does not match info/block_name")
        for c1, b1, c2, b2 in self.block2_pairs:
            if not (0 <= c1 < n_c and 0 <= c2 < n_c and 0 <= b1 < len(self.G) and 0 <= b2 < len(self.G)):
                raise X0SourceError("x0_source/block2_pairs has an out-of-range index")
        if len(self.inner2_pairs) != n_inner:
            raise X0SourceError("x0_source/inner2_pairs does not match info/inner_name")
        if np.any(self.inner2_pairs < 0) or np.any(self.inner2_pairs >= n_in):
            raise X0SourceError("x0_source/inner2_pairs has an out-of-range index")
        if not np.isclose(beta, info_beta, rtol=1e-12, atol=0.0):
            raise X0SourceError("x0_source beta does not match info/beta")
        if self.iw_offset != self.n_iwG // 2:
            raise X0SourceError("x0_source iw_offset is inconsistent with n_iwG")
        for wb in self.list_wb:
            if wb < 0 or self.iw_offset - self.num_wf < 0 or self.iw_offset + self.num_wf + wb > self.n_iwG:
                raise X0SourceError(f"bosonic index {wb} is outside the stored G frequency window")
            if f"w{wb}" not in x0_loc_wb:
                raise X0SourceError(f"input/X0_loc is missing for bosonic index {wb}")

    def check_target(self, wb, q_label):
        """Raise X0SourceError if (wb, q_label) cannot be reconstructed."""
        if int(wb) not in self.list_wb:
            raise X0SourceError(f"bosonic index {wb} is not stored (available: {self.list_wb})")
        parse_q_label(q_label, self.div)

    def get(self, wb, q_label):
        """X0_q in h5bse format: {(b_l, b_r): array [n_in2, n_in2, 2*num_wf]}."""
        wb = int(wb)
        self.check_target(wb, q_label)
        q = parse_q_label(q_label, self.div)
        if q not in self._kq_cache:
            if len(self._kq_cache) >= _KQ_CACHE_SIZE:
                self._kq_cache.pop(next(iter(self._kq_cache)))
            self._kq_cache[q] = k_minus_q_index(self.div, q)
        return x0_from_g(self.G, self.block2_pairs, self.inner2_pairs, wb, self.num_wf,
                         self.iw_offset, self._kq_cache[q])


class X0SourceWriter:
    """Write the ``x0_source`` group (rank 0 only).

    Usage: construct (removes stale X0 data), write_meta, write_g_slice for every
    slice, then finish(). Gbar (k-averaged G) is accumulated for X0_loc.
    """

    def __init__(self, h5f, groupname="", compression=None, compression_opts=None):
        root = bse_root(groupname)
        for path in (root + "/input/X0_q", root + "/" + GROUP):
            if path in h5f:
                del h5f[path]
        self._grp = h5f.create_group(root + "/" + GROUP)
        self._grp.attrs["format"] = FORMAT
        self._grp.attrs["version"] = VERSION
        self._grp.attrs["complete"] = 0
        self._compression = compression
        self._compression_opts = compression_opts
        self.Gbar = None

    def write_meta(self, block1_names, block2_pairs, inner2_pairs, num_wf, list_wb, n_iwG, div, beta,
                   shape_c_in):
        """shape_c_in = (n_c, N_in)."""
        g = self._grp
        g.attrs["num_wf"] = int(num_wf)
        g.attrs["list_wb"] = np.asarray(list_wb, dtype=int)
        g.attrs["n_iwG"] = int(n_iwG)
        g.attrs["iw_offset"] = int(n_iwG) // 2
        g.attrs["div"] = np.asarray(div, dtype=int)
        g.attrs["beta"] = float(beta)
        g.create_dataset("block1_names", data=np.array([str(n) for n in block1_names], dtype=object),
                         dtype=h5py.string_dtype())
        g.create_dataset("block2_pairs", data=np.asarray(block2_pairs, dtype=np.int64).reshape(-1, 4))
        g.create_dataset("inner2_pairs", data=np.asarray(inner2_pairs, dtype=np.int64).reshape(-1, 2))
        n_c, n_in = shape_c_in
        n_k = int(np.prod(div))
        shape = (n_c, n_c, n_in, n_in, int(n_iwG), n_k)
        kwargs = {}
        if self._compression is not None:
            kwargs["compression"] = self._compression
            if self._compression == "gzip" and self._compression_opts is not None:
                kwargs["compression_opts"] = self._compression_opts
        self._dsets = [
            g.create_dataset(f"G/{b}", shape=shape, dtype=np.complex128,
                             chunks=_g_chunks(int(n_iwG), n_k), **kwargs)
            for b in range(len(block1_names))
        ]
        self.Gbar = [np.zeros(shape[:5], dtype=np.complex128) for _ in block1_names]

    def write_g_slice(self, b, c1, c2, i1, i2, g_wk):
        """g_wk: array [iw, k] for single-particle block index b."""
        self._dsets[b][c1, c2, i1, i2, :, :] = g_wk
        self.Gbar[b][c1, c2, i1, i2, :] = np.mean(g_wk, axis=1)

    def finish(self):
        self._grp.attrs["complete"] = 1
