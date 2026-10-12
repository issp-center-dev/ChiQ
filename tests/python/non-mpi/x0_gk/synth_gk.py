"""Synthetic data and independent reference implementations for x0_gk tests."""

from itertools import product

import h5py
import numpy as np

from chiq import h5bse
from chiq.x0_gk import X0SourceWriter


def make_pairs(n_c, n_b1, n_in, nonlocal_pairs=True):
    block2 = [(c1, b1, c2, b2)
              for c1, c2, b1, b2 in product(range(n_c), range(n_c), range(n_b1), range(n_b1))
              if nonlocal_pairs or c1 == c2]
    inner2 = [(i1, i2) for i1, i2 in product(range(n_in), repeat=2)]
    return block2, inner2


def random_g(rng, n_b1, n_c, n_in, n_iwG, n_k):
    shape = (n_c, n_c, n_in, n_in, n_iwG, n_k)
    return [rng.normal(size=shape) + 1j * rng.normal(size=shape) for _ in range(n_b1)]


def model_g(beta, n_b1, n_c, n_in, n_iwG, div, seed=0):
    """Non-interacting G(iw, k) = (iw - H(k))^-1 for a random tight-binding H(k)
    without inversion symmetry (so that X0(q) != X0(-q))."""
    rng = np.random.default_rng(seed)
    dim = n_c * n_in
    n_k = int(np.prod(div))
    iw = 1j * np.pi * (2 * (np.arange(n_iwG) - n_iwG // 2) + 1) / beta
    kvec = np.indices(div).reshape(3, -1).T * (2 * np.pi / np.asarray(div))
    G = []
    for _ in range(n_b1):
        h0 = rng.normal(size=(dim, dim)) + 1j * rng.normal(size=(dim, dim))
        h0 = 0.5 * (h0 + h0.conj().T)
        t = [rng.normal(size=(dim, dim)) + 1j * rng.normal(size=(dim, dim)) for _ in range(3)]
        g = np.empty((n_c, n_c, n_in, n_in, n_iwG, n_k), dtype=complex)
        for ik, k in enumerate(kvec):
            hk = h0.copy()
            for l in range(3):
                hop = t[l] * np.exp(1j * k[l])
                hk += 0.3 * (hop + hop.conj().T)
            for iiw, w in enumerate(iw):
                gk = np.linalg.inv(w * np.eye(dim) - hk).reshape(n_c, n_in, n_c, n_in)
                g[:, :, :, :, iiw, ik] = gk.transpose(0, 2, 1, 3)
        G.append(g)
    return G


def flat_k(k, div):
    return k[2] + div[2] * (k[1] + div[1] * k[0])


def oracle_x0(G, block2, inner2, wb, num_wf, iw_offset, div, q, sign=-1):
    """Explicit loops: X0 = -(1/n_k) sum_k G_l(t, k) G_r(t + wb, k + sign*q)."""
    n_k = int(np.prod(div))
    lo = iw_offset - num_wf
    kq = []
    for k in product(range(div[0]), range(div[1]), range(div[2])):
        kq.append(flat_k(tuple((k[l] + sign * q[l]) % div[l] for l in range(3)), div))
    out = {}
    for b_l, (c1, b1, c2, b2) in enumerate(block2):
        for b_r, (c3, b3, c4, b4) in enumerate(block2):
            if b1 != b3 or b2 != b4:
                continue
            x = np.zeros((len(inner2), len(inner2), 2 * num_wf), dtype=complex)
            for i, (i1, i2) in enumerate(inner2):
                for j, (i3, i4) in enumerate(inner2):
                    for t in range(2 * num_wf):
                        s = 0j
                        for k in range(n_k):
                            s += G[b1][c3, c1, i3, i1, lo + t, k] * G[b2][c2, c4, i2, i4, lo + t + wb, kq[k]]
                        x[i, j, t] = -s / n_k
            out[(b_l, b_r)] = x
    return out


def fft_x0_all_q(G, block2, inner2, wb, num_wf, iw_offset, div):
    """Same formula as SumkDFTChi.calc_chi0_fft, for all q. Returns {pair: [i, j, t, q]}."""
    n_k = int(np.prod(div))
    lo = iw_offset - num_wf
    n_wf2 = 2 * num_wf
    minus_r = np.array([flat_k(tuple((-r[l]) % div[l] for l in range(3)), div)
                        for r in product(range(div[0]), range(div[1]), range(div[2]))])

    def k2r(gwk):
        return np.fft.ifftn(gwk.reshape((gwk.shape[0],) + tuple(div)), axes=(1, 2, 3)).reshape(gwk.shape[0], n_k)

    def r2k(gwr):
        return np.fft.fftn(gwr.reshape((gwr.shape[0],) + tuple(div)), axes=(1, 2, 3)).reshape(gwr.shape[0], n_k)

    out = {}
    for b_l, (c1, b1, c2, b2) in enumerate(block2):
        for b_r, (c3, b3, c4, b4) in enumerate(block2):
            if b1 != b3 or b2 != b4:
                continue
            x = np.zeros((len(inner2), len(inner2), n_wf2, n_k), dtype=complex)
            for i, (i1, i2) in enumerate(inner2):
                for j, (i3, i4) in enumerate(inner2):
                    g1 = G[b1][c3, c1, i3, i1, lo:lo + n_wf2 + wb, :]
                    g2 = G[b2][c2, c4, i2, i4, lo:lo + n_wf2 + wb, :]
                    g1r, g2r = k2r(g1), k2r(g2)
                    chi_r = -g1r[0:n_wf2, :] * g2r[:, minus_r][wb:wb + n_wf2, :]
                    x[i, j] = r2k(chi_r)
            out[(b_l, b_r)] = x
    return out


def write_info(path, block_names, inner_names, beta, groupname=""):
    bs = h5bse.h5BSE(path, groupname)
    bs.save(key=("block_name",), data=list(block_names))
    bs.save(key=("inner_name",), data=list(inner_names))
    bs.save(key=("beta",), data=beta)


def write_gk_file(path, G, block2, inner2, num_wf, list_wb, div, beta, groupname=""):
    """A gk-mode input file: info, x0_source, X0_loc (from Gbar)."""
    from chiq.x0_gk import x0_loc_from_gbar

    n_b1 = len(G)
    n_c, n_in, n_iwG = G[0].shape[0], G[0].shape[2], G[0].shape[4]
    write_info(path, [str(p) for p in block2], [str(p) for p in inner2], beta, groupname)
    with h5py.File(path, "a") as h5f:
        w = X0SourceWriter(h5f, groupname, compression="gzip", compression_opts=4)
        w.write_meta([f"b{b}" for b in range(n_b1)], block2, inner2, num_wf, list_wb, n_iwG, div, beta,
                     (n_c, n_in))
        for b in range(n_b1):
            for c1, c2, i1, i2 in product(range(n_c), range(n_c), range(n_in), range(n_in)):
                w.write_g_slice(b, c1, c2, i1, i2, G[b][c1, c2, i1, i2])
        gbar = w.Gbar
        w.finish()
    bs = h5bse.h5BSE(path, groupname)
    for wb in list_wb:
        bs.save(key=("X0_loc", wb), data=x0_loc_from_gbar(gbar, block2, inner2, wb, num_wf, n_iwG // 2))
