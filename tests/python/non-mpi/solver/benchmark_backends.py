#!/usr/bin/env python3
"""Speed comparison of the ChiQ solver backends (cpp vs numpy).

Not a pytest test (no ``test_`` prefix) -- run it directly:

    PYTHONPATH="$REPO/python/package:$REPO/build/src" \
        python3 tests/python/non-mpi/solver/benchmark_backends.py

For each calc type it builds well-conditioned random inputs of a representative
size ONCE and times both backends on the *same* inputs, reporting the per-call
wall time of each backend and the numpy/cpp ratio. The cpp backend uses the
compiled (optionally OpenMP-threaded) C++ solver; the numpy backend is a
single-threaded pure-Python reimplementation, so it is expected to be slower --
this script quantifies by how much.
"""

import sys
import time
import numpy as np

from chiq.solver import get_solver

RNG = np.random.default_rng(0)


def _info(nb, nw, n_in):
    return ([n_in * nw] * nb, [n_in] * (nb * nw), [n_in] * nb)


def _rand_c(shape):
    return RNG.standard_normal(shape) + 1j * RNG.standard_normal(shape)


def _wc(n):
    """A strictly diagonally-dominant (hence invertible) complex matrix.

    Setting each diagonal entry's magnitude above that row's off-diagonal
    absolute sum guarantees invertibility for any n / seed (Levy-Desplanques).
    """
    m = _rand_c((n, n))
    off = np.abs(m).sum(axis=1) - np.abs(np.diag(m))
    np.fill_diagonal(m, off + 1.0)
    return m


def _c_full(nb, n_in):
    return {(i, j): _rand_c((n_in, n_in)) for i in range(nb) for j in range(nb)}


def _c_wc(nb, n_in):
    return {(i, i): _wc(n_in) for i in range(nb)}


def _b_wc(nb, nw, n_in):
    # B-type: each same-iw super-block is block-diagonally dominant (invertible)
    # -- strong diagonal blocks, small off-diagonal blocks.
    out = {}
    for iw in range(nw):
        for b in range(nb):
            for b2 in range(nb):
                if b == b2:
                    out[(iw * nb + b, iw * nb + b2)] = _wc(n_in) + n_in * np.eye(n_in)
                else:
                    out[(iw * nb + b, iw * nb + b2)] = 0.1 * _rand_c((n_in, n_in))
    return out


def _b_full(nb, nw, n_in):
    return {(iw * nb + b, iw * nb + b2): _rand_c((n_in, n_in))
            for iw in range(nw) for b in range(nb) for b2 in range(nb)}


def _a_full(nb, nw, n_in):
    return {(i, j): _rand_c((n_in * nw, n_in * nw))
            for i in range(nb) for j in range(nb)}


def _inputs(calc, nb, nw, n_in):
    if calc == "chi0":
        b = _b_full(nb, nw, n_in)
        return {"chi0_loc": _c_full(nb, n_in), "X0_loc": b,
                "X0_q": {k: _rand_c((n_in, n_in)) for k in b}}
    if calc == "rpa":
        return {"chi0_q": _c_full(nb, n_in),
                "gamma0": {(i, i): 0.01 * _rand_c((n_in, n_in)) for i in range(nb)}}
    if calc == "rrpa":
        return {"chi0_loc": _c_wc(nb, n_in), "chi_loc": _c_wc(nb, n_in),
                "chi0_q": _c_full(nb, n_in)}
    if calc == "bse":
        return {"X0_loc": _b_wc(nb, nw, n_in), "X0_q": _b_wc(nb, nw, n_in),
                "X_loc": _a_full(nb, nw, n_in), "chi_loc": _c_wc(nb, n_in)}
    if calc == "scl":
        return {"X0_loc": _b_wc(nb, nw, n_in), "X0_q": _b_wc(nb, nw, n_in),
                "Phi": _b_full(nb, nw, n_in), "Phi_sum": _c_wc(nb, n_in)}
    raise ValueError(calc)


def _time_calc(backend, calc, inputs, info, repeat):
    """Best-of-`repeat` wall time of one calc() on the *given* inputs.

    A warm-up call precedes the timed runs so first-use initialization is not
    charged to either backend. set() is excluded from the timed region.
    (set() copies its inputs, so the same `inputs` dict is reused for every
    backend/run without cross-contamination.)
    """
    iA, iB, iC = info
    best = float("inf")
    for r in range(repeat + 1):  # r == 0 is an untimed warm-up
        s = get_solver(backend, 12.0, iA, iB, iC)
        for name, bm in inputs.items():
            s.set(bm, name)
        t0 = time.perf_counter()
        s.calc(calc)
        dt = time.perf_counter() - t0
        if r > 0:
            best = min(best, dt)
    return best


def main():
    # (nb, nw, n_in): a couple of representative single-q solves.
    sizes = [(1, 20, 2), (2, 20, 2), (2, 40, 2)]
    calcs = ["chi0", "rpa", "rrpa", "bse", "scl"]
    repeat = 3
    print(f"per-call best-of-{repeat} wall time (one (omega,q) solve), seconds")
    print("both backends timed on identical inputs; warm-up excluded\n")
    header = f"{'calc':6} {'nb,nw,n_in':>12}  {'cpp [s]':>10} {'numpy [s]':>10} {'numpy/cpp':>10}"
    print(header)
    print("-" * len(header))
    for nb, nw, n_in in sizes:
        info = _info(nb, nw, n_in)
        for calc in calcs:
            inputs = _inputs(calc, nb, nw, n_in)  # built ONCE, shared by both backends
            t_cpp = _time_calc("cpp", calc, inputs, info, repeat)
            t_np = _time_calc("numpy", calc, inputs, info, repeat)
            ratio = t_np / t_cpp if t_cpp > 0 else float("inf")
            print(f"{calc:6} {f'{nb},{nw},{n_in}':>12}  "
                  f"{t_cpp:10.4g} {t_np:10.4g} {ratio:9.1f}x")
        print()
    print("Note: cpp is the production backend; numpy is single-threaded and "
          "reproduces the same results (see test_backend_agreement.py).")


if __name__ == "__main__":
    sys.exit(main())
