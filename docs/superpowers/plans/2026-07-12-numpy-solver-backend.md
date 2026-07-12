# NumPy Solver Backend Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a `backend = "cpp" | "numpy" | "cupy"` switch to the ChiQ BSE solver, with a NumPy backend that numerically matches the existing C++ backend, so later phases can add a CuPy/GPU backend and IR basis without touching driver logic.

**Architecture:** A new `chiq.solver` subpackage with three separated concerns — `layout` (block-dict ↔ dense block-matrix algebra, the sole owner of the A/B/C index conventions), `kernels` (the five susceptibility calculations written against an injected array module `xp`), and a stateful `set`/`calc`/`get` facade (`CppSolver` wrapping the existing pybind module, `NumpySolver` running the kernels, `CupySolver` skeleton). `chiq_main.py` obtains its solver only through a factory.

**Tech Stack:** Python ≥3.10, NumPy, the existing pybind11 module `bse_solver`, pytest. The block-matrix data model is `dict[(int,int) -> np.ndarray(complex128)]`.

## Global Constraints

- **Zero behavior change by default.** Default backend is `cpp`; existing runs/outputs are unchanged. `chiq_main.py` currently constructs `bse_solver.BSESolver` inline — after Task 7 it goes through `chiq.solver.get_solver`, and the `cpp` path must be byte-for-byte equivalent.
- **Compatibility target = mathematical equivalence within tolerance**, not bitwise. Agreement metric (Task 12): elementwise `np.testing.assert_allclose(numpy, cpp, rtol=1e-8, atol=1e-10, equal_nan=True)` per block, preceded by exact NaN and signed-infinity position masks.
- **The NumPy kernels mirror the C++ operation sequence verbatim** (explicit `block_inverse` then `matmul` wherever `bse.hpp` writes `inverse(...) * ...`), to agree with the C++ oracle at the tolerance above. The `solve`-instead-of-`inverse` optimization from the spec §3.5 is deferred to the later CuPy/performance phase.
- **Uniform-dimension invariant:** every block has the same inner dimension `n_in` and frequency count `nw` (`chiq_main.py:524-526`; `sumk_dft_chi.py:579-580`). The layout module requires this and raises on non-uniform matrix-info. Heterogeneous blocks are unsupported.
- **dtype:** inputs normalized to `complex128` on `set`; `get` returns `complex128`.
- **Backend names accepted by `calc`** are case-tolerant: `"BSE"|"bse"`, `"RPA"|"rpa"`, `"RRPA"|"rrpa"`, `"SCL"|"scl"`, `"chi0"`.
- **`get` before compute returns `{}`** (empty dict), never an exception. Unknown set/get/calc names raise `RuntimeError("Invalid type '<name>'")`.
- **Missing-input error strings** use the C++ *internal* names — missing `X0_loc` reports `'X0_Loc'` (capital L); others match their set names.
- **`I_q` and `I_q_scl` are one result slot** (both map to the single `Iq` member).
- **Spec:** `docs/superpowers/specs/2026-07-12-numpy-backend-and-packaging-design.md`. Packaging (§5, scikit-build-core, `_bse_solver` rename, entry points) is a **separate plan** — do not do it here; keep importing the current top-level `bse_solver` module.

## File Structure

- Create `python/package/chiq/solver/__init__.py` — `get_solver` factory, re-exports.
- Create `python/package/chiq/solver/base.py` — `SolverBase` abstract facade + shared name/type tables.
- Create `python/package/chiq/solver/layout.py` — block-matrix algebra: `parse_matrix_info`, `add/sub/scale/identity`, `connected_components`, `block_inverse`, `matmul`, `sumfreq_a`, `sumfreq_b`, `convert_a2b`, `convert_b2a`, `prod_ab`.
- Create `python/package/chiq/solver/kernels.py` — `calc_chi0`, `calc_bse`, `calc_rpa`, `calc_rrpa`, `calc_scl`, each `(state, beta, dims) -> dict[name -> bm]`.
- Create `python/package/chiq/solver/cpp.py` — `CppSolver`.
- Create `python/package/chiq/solver/numpy.py` — `NumpySolver`.
- Create `python/package/chiq/solver/cupy.py` — `CupySolver` skeleton.
- Modify `python/scripts/chiq_main.py` — construct the solver via `get_solver`.
- Modify `python/package/chiq/bse_toml.py` — parse+validate `[chiq_main] backend`.
- Create tests under `tests/python/non-mpi/solver/` — unit tests for layout/kernels/facade.
- Modify existing fixtures `tests/python/non-mpi/bsetool_{BSE,RPA,SCL,RRPA}/test_*.py` — parametrize over backend (Task 12).

Layout conventions (all reproduced from `src/bse.hpp`, `src/block_matrix.hpp`):
`nb` = number of C-vertices, `nw` = number of frequencies, `n_in` = inner dim.
- C-type: key `b∈[0,nb)`, block `n_in×n_in`, dims `[n_in]*nb`.
- B-type: key `iw*nb+b`, block `n_in×n_in`, dims `[n_in]*(nb*nw)`. Inputs are iw-diagonal; intermediates (from `convert_a2b`) may have cross-iw keys.
- A-type: key `b∈[0,nb)`, block `(n_in*nw)×(n_in*nw)`, inner flattened `in*nw+iw`, dims `[n_in*nw]*nb`.
- Missing key = structural zero.

---

### Task 1: layout — matrix-info parsing and validation

**Files:**
- Create: `python/package/chiq/solver/__init__.py` (empty for now: `# chiq solver subpackage`)
- Create: `python/package/chiq/solver/layout.py`
- Test: `tests/python/non-mpi/solver/test_layout_info.py`

**Interfaces:**
- Produces: `parse_matrix_info(info_A, info_B, info_C) -> (nb, nw, n_in)`. Raises `ValueError` on non-uniform or inconsistent info. `info_*` are sequences of ints.

- [ ] **Step 1: Write the failing test**

```python
# tests/python/non-mpi/solver/test_layout_info.py
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
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python -m pytest tests/python/non-mpi/solver/test_layout_info.py -v`
Expected: FAIL (module `chiq.solver.layout` has no `parse_matrix_info`).

- [ ] **Step 3: Write minimal implementation**

```python
# python/package/chiq/solver/layout.py
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
```

- [ ] **Step 4: Run test to verify it passes**

Run: `python -m pytest tests/python/non-mpi/solver/test_layout_info.py -v`
Expected: PASS (4 tests). Note: to import `chiq`, run from a tree where `python/package` is importable — source the build `chiqvars.sh`, or `PYTHONPATH=python/package python -m pytest ...`.

- [ ] **Step 5: Commit**

```bash
git add python/package/chiq/solver/__init__.py python/package/chiq/solver/layout.py tests/python/non-mpi/solver/test_layout_info.py
git commit -m "feat(solver): layout.parse_matrix_info with uniform-dim validation"
```

---

### Task 2: layout — elementwise add/sub/scale and identity

**Files:**
- Modify: `python/package/chiq/solver/layout.py`
- Test: `tests/python/non-mpi/solver/test_layout_elementwise.py`

**Interfaces:**
- Consumes: nothing.
- Produces: `add(a, b) -> dict`, `sub(a, b) -> dict`, `scale(s, bm) -> dict`, `identity(dims) -> dict`. Block matrices are `dict[(int,int)->ndarray]`; `dims` is a sequence of per-vertex dimensions.

- [ ] **Step 1: Write the failing test**

```python
# tests/python/non-mpi/solver/test_layout_elementwise.py
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
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python -m pytest tests/python/non-mpi/solver/test_layout_elementwise.py -v`
Expected: FAIL (`add`/`sub`/`scale`/`identity` undefined).

- [ ] **Step 3: Write minimal implementation**

```python
# append to python/package/chiq/solver/layout.py

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
```

- [ ] **Step 4: Run test to verify it passes**

Run: `python -m pytest tests/python/non-mpi/solver/test_layout_elementwise.py -v`
Expected: PASS (3 tests).

- [ ] **Step 5: Commit**

```bash
git add python/package/chiq/solver/layout.py tests/python/non-mpi/solver/test_layout_elementwise.py
git commit -m "feat(solver): layout add/sub/scale/identity"
```

---

### Task 3: layout — connected components and block inverse (with fill-in)

**Files:**
- Modify: `python/package/chiq/solver/layout.py`
- Test: `tests/python/non-mpi/solver/test_layout_inverse.py`

**Interfaces:**
- Consumes: nothing.
- Produces: `connected_components(bm, nvert) -> list[list[int]]` (components under edge `(i,j) or (j,i) present`, each component's vertices sorted, components sorted by first vertex). `block_inverse(bm, dims) -> dict` (per-component dense inverse; emits **every** intra-component vertex pair, numerical zeros retained; cross-component pairs absent).

- [ ] **Step 1: Write the failing test**

```python
# tests/python/non-mpi/solver/test_layout_inverse.py
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
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python -m pytest tests/python/non-mpi/solver/test_layout_inverse.py -v`
Expected: FAIL (`connected_components`/`block_inverse` undefined).

- [ ] **Step 3: Write minimal implementation**

```python
# append to python/package/chiq/solver/layout.py

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
```

- [ ] **Step 4: Run test to verify it passes**

Run: `python -m pytest tests/python/non-mpi/solver/test_layout_inverse.py -v`
Expected: PASS (3 tests).

- [ ] **Step 5: Commit**

```bash
git add python/package/chiq/solver/layout.py tests/python/non-mpi/solver/test_layout_inverse.py
git commit -m "feat(solver): layout connected_components + block_inverse with fill-in"
```

---

### Task 4: layout — block matrix product

**Files:**
- Modify: `python/package/chiq/solver/layout.py`
- Test: `tests/python/non-mpi/solver/test_layout_matmul.py`

**Interfaces:**
- Consumes: nothing.
- Produces: `matmul(a, b, dims) -> dict`. Result key `(i,j)` present **iff ∃k with `(i,k)` in `a` and `(k,j)` in `b`** (structural, value-independent, matching `block_matrix.hpp operator*`); value `Σ_k a[(i,k)] @ b[(k,j)]`. `dims` gives `nvert=len(dims)`.

- [ ] **Step 1: Write the failing test**

```python
# tests/python/non-mpi/solver/test_layout_matmul.py
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
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python -m pytest tests/python/non-mpi/solver/test_layout_matmul.py -v`
Expected: FAIL (`matmul` undefined).

- [ ] **Step 3: Write minimal implementation**

```python
# append to python/package/chiq/solver/layout.py

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
```

- [ ] **Step 4: Run test to verify it passes**

Run: `python -m pytest tests/python/non-mpi/solver/test_layout_matmul.py -v`
Expected: PASS (2 tests).

- [ ] **Step 5: Commit**

```bash
git add python/package/chiq/solver/layout.py tests/python/non-mpi/solver/test_layout_matmul.py
git commit -m "feat(solver): layout.matmul with structural key rule"
```

---

### Task 5: layout — Sumfreq_A and Sumfreq_B

**Files:**
- Modify: `python/package/chiq/solver/layout.py`
- Test: `tests/python/non-mpi/solver/test_layout_sumfreq.py`

**Interfaces:**
- Consumes: nothing.
- Produces: `sumfreq_a(bm, nb, nw, n_in) -> dict` (A→C; only keys present in `bm`; `out[(i,j)][a,c] = Σ_{iw1,iw2} A[(i,j)][a*nw+iw1, c*nw+iw2]`). `sumfreq_b(bm, nb, nw, n_in) -> dict` (B→C; **all** `nb²` keys emitted; `out[(i,j)] = Σ_iw B[(iw*nb+i, iw*nb+j)]` over present keys).

- [ ] **Step 1: Write the failing test**

```python
# tests/python/non-mpi/solver/test_layout_sumfreq.py
import numpy as np
from chiq.solver import layout

def test_sumfreq_b_sums_over_frequency_and_emits_all_keys():
    nb, nw, n_in = 2, 2, 1
    bm = {}
    # only diagonal (i==j) same-iw keys present
    bm[(0 * nb + 0, 0 * nb + 0)] = np.array([[1.0]], dtype=complex)
    bm[(1 * nb + 0, 1 * nb + 0)] = np.array([[3.0]], dtype=complex)  # iw=1, b=0
    out = layout.sumfreq_b(bm, nb, nw, n_in)
    assert set(out) == {(i, j) for i in range(nb) for j in range(nb)}
    assert np.allclose(out[(0, 0)], [[4.0]])  # 1 + 3
    assert np.allclose(out[(0, 1)], [[0.0]])

def test_sumfreq_a_sums_over_both_frequency_indices():
    nb, nw, n_in = 1, 2, 1
    # A block is (n_in*nw)=(2) x 2; inner flattened as in*nw+iw -> here in=0 so index==iw
    A = np.array([[1, 2], [3, 4]], dtype=complex)
    bm = {(0, 0): A}
    out = layout.sumfreq_a(bm, nb, nw, n_in)
    assert set(out) == {(0, 0)}
    assert np.allclose(out[(0, 0)], [[1 + 2 + 3 + 4]])
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python -m pytest tests/python/non-mpi/solver/test_layout_sumfreq.py -v`
Expected: FAIL (`sumfreq_a`/`sumfreq_b` undefined).

- [ ] **Step 3: Write minimal implementation**

```python
# append to python/package/chiq/solver/layout.py

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
```

- [ ] **Step 4: Run test to verify it passes**

Run: `python -m pytest tests/python/non-mpi/solver/test_layout_sumfreq.py -v`
Expected: PASS (2 tests).

- [ ] **Step 5: Commit**

```bash
git add python/package/chiq/solver/layout.py tests/python/non-mpi/solver/test_layout_sumfreq.py
git commit -m "feat(solver): layout sumfreq_a / sumfreq_b"
```

---

### Task 6: layout — Convert_A2B, Convert_B2A, prod_ab

**Files:**
- Modify: `python/package/chiq/solver/layout.py`
- Test: `tests/python/non-mpi/solver/test_layout_convert.py`

**Interfaces:**
- Consumes: `matmul`.
- Produces: `convert_a2b(bm, nb, nw, n_in) -> dict` (A→B; key `(iw1*nb+b1, iw2*nb+b2)` emitted iff the extracted block is not exactly zero — `isZero(0)` pruning). `convert_b2a(bm, nb, nw, n_in) -> dict` (B→A; key `(b1,b2)` emitted iff its assembled A block is not exactly zero). `prod_ab(a_mat, b_mat, nb, nw, n_in) -> dict` = `convert_b2a(matmul(convert_a2b(a_mat), b_mat, [n_in]*(nb*nw)))` (A-type result), reproducing `bse.hpp Prod_A_B`.

- [ ] **Step 1: Write the failing test**

```python
# tests/python/non-mpi/solver/test_layout_convert.py
import numpy as np
from chiq.solver import layout

def test_a2b_b2a_roundtrip():
    nb, nw, n_in = 1, 2, 1
    A = np.array([[1, 2], [3, 4]], dtype=complex)  # in*nw+iw indexing, n_in=1
    a = {(0, 0): A}
    b = layout.convert_a2b(a, nb, nw, n_in)
    # keys are (iw1*nb+0, iw2*nb+0) for nonzero extracted entries
    assert (0, 0) in b and (1, 1) in b  # A[0,0]=1, A[1,1]=4
    back = layout.convert_b2a(b, nb, nw, n_in)
    assert np.allclose(back[(0, 0)], A)

def test_a2b_prunes_exact_zero_blocks():
    nb, nw, n_in = 1, 2, 1
    A = np.array([[1, 0], [0, 4]], dtype=complex)
    b = layout.convert_a2b({(0, 0): A}, nb, nw, n_in)
    assert set(b) == {(0, 0), (1, 1)}  # cross-freq (0,1),(1,0) are exactly zero -> pruned

def test_prod_ab_matches_dense():
    nb, nw, n_in = 1, 2, 1
    A = np.array([[1, 2], [3, 4]], dtype=complex)
    # B-type diagonal in iw: keys (iw,iw)
    Bmat = {(0, 0): np.array([[2]], dtype=complex), (1, 1): np.array([[5]], dtype=complex)}
    got = layout.prod_ab(A_as := {(0, 0): A}, Bmat, nb, nw, n_in)
    # dense: convert A to B (full 2x2 in freq), multiply by diag(2,5), convert back
    Bdense = np.diag([2.0, 5.0]).astype(complex)
    assert np.allclose(got[(0, 0)], A @ Bdense)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python -m pytest tests/python/non-mpi/solver/test_layout_convert.py -v`
Expected: FAIL (`convert_a2b`/`convert_b2a`/`prod_ab` undefined).

- [ ] **Step 3: Write minimal implementation**

```python
# append to python/package/chiq/solver/layout.py

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
    exactly zero (matches Convert_B2A isZero(0) pruning — gate on the
    assembled VALUE, not on input key presence: matmul can produce a
    present-but-exactly-zero B block via cancellation)."""
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
```

- [ ] **Step 4: Run test to verify it passes**

Run: `python -m pytest tests/python/non-mpi/solver/test_layout_convert.py -v`
Expected: PASS (3 tests).

- [ ] **Step 5: Commit**

```bash
git add python/package/chiq/solver/layout.py tests/python/non-mpi/solver/test_layout_convert.py
git commit -m "feat(solver): layout convert_a2b/convert_b2a/prod_ab"
```

---

### Task 7: SolverBase + factory + CppSolver; route chiq_main and bse_toml

**Files:**
- Create: `python/package/chiq/solver/base.py`
- Create: `python/package/chiq/solver/cpp.py`
- Modify: `python/package/chiq/solver/__init__.py`
- Modify: `python/package/chiq/bse_toml.py`
- Modify: `python/scripts/chiq_main.py:551` (solver construction) and its solver calls
- Test: `tests/python/non-mpi/solver/test_factory_cpp.py`

**Interfaces:**
- Consumes: existing `bse_solver.BSESolver`.
- Produces: `chiq.solver.get_solver(backend, beta, info_A, info_B, info_C) -> SolverBase`; `SolverBase` with `set(bm, name)`, `calc(calc_type)`, `get(name) -> dict`. `bse_toml.load_params_from_toml` returns a `dict_tool` that now includes `"backend"` (default `"cpp"`).

- [ ] **Step 1: Write the failing test**

```python
# tests/python/non-mpi/solver/test_factory_cpp.py
import numpy as np
import pytest
from chiq.solver import get_solver

def _info(nb, nw, n_in):
    return ([n_in * nw] * nb, [n_in] * (nb * nw), [n_in] * nb)

def test_factory_returns_cpp_by_default_and_runs_chi0():
    nb, nw, n_in = 1, 2, 1
    iA, iB, iC = _info(nb, nw, n_in)
    solver = get_solver("cpp", beta=10.0, info_A=iA, info_B=iB, info_C=iC)
    # minimal chi0 inputs
    solver.set({(0, 0): np.zeros((n_in, n_in), dtype=complex)}, "chi0_loc")
    solver.set({(0, 0): np.zeros((n_in, n_in), dtype=complex), (1, 1): np.zeros((n_in, n_in), dtype=complex)}, "X0_loc")
    solver.set({(0, 0): np.zeros((n_in, n_in), dtype=complex), (1, 1): np.zeros((n_in, n_in), dtype=complex)}, "X0_q")
    solver.calc("chi0")
    out = solver.get("chi0_q")
    assert isinstance(out, dict)

def test_factory_rejects_unknown_backend():
    iA, iB, iC = _info(1, 2, 1)
    with pytest.raises(ValueError):
        get_solver("gpu999", 10.0, iA, iB, iC)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python -m pytest tests/python/non-mpi/solver/test_factory_cpp.py -v` (needs the built `bse_solver` importable — source `chiqvars.sh` from the build dir first).
Expected: FAIL (`get_solver` undefined).

- [ ] **Step 3: Write minimal implementation**

```python
# python/package/chiq/solver/base.py
"""Solver facade base class and shared name tables."""

import abc

# set-name -> layout type
SET_LAYOUT = {
    "X_loc": "A", "X0_loc": "B", "X0_q": "B", "Phi": "B",
    "chi_loc": "C", "chi0_loc": "C", "chi0_q": "C", "gamma0": "C", "Phi_sum": "C",
}
CALC_TYPES = {"chi0", "bse", "rpa", "rrpa", "scl"}
# valid get() output names (from bse_solver_pybind.cpp getMatrix)
GET_NAMES = {"chi0_q", "chi_q", "chi_q_rpa", "chi_q_rrpa", "chi_q_scl", "I_q", "I_q_scl"}


class SolverBase(abc.ABC):
    def __init__(self, beta, info_A, info_B, info_C):
        self.beta = float(beta)
        self.info_A = list(int(x) for x in info_A)
        self.info_B = list(int(x) for x in info_B)
        self.info_C = list(int(x) for x in info_C)

    @abc.abstractmethod
    def set(self, bm, name):
        ...

    @abc.abstractmethod
    def calc(self, calc_type):
        ...

    @abc.abstractmethod
    def get(self, name):
        ...
```

```python
# python/package/chiq/solver/cpp.py
"""C++ backend: thin wrapper over the pybind module `bse_solver`."""

import numpy as np
import bse_solver  # top-level pybind module (packaging phase renames to chiq._bse_solver)

from .base import SolverBase


class CppSolver(SolverBase):
    def __init__(self, beta, info_A, info_B, info_C):
        super().__init__(beta, info_A, info_B, info_C)
        self._impl = bse_solver.BSESolver(
            self.beta,
            np.array(self.info_A),
            np.array(self.info_B),
            np.array(self.info_C),
        )

    def set(self, bm, name):
        self._impl.set(bm, name)

    def calc(self, calc_type):
        self._impl.calc(calc_type)

    def get(self, name):
        return self._impl.get(name)
```

```python
# python/package/chiq/solver/__init__.py
"""chiq solver subpackage: backend factory."""

from .base import SolverBase


def get_solver(backend, beta, info_A, info_B, info_C):
    """Return a solver facade for the named backend.

    backend: "cpp" (default in the driver), "numpy", or "cupy".
    """
    name = str(backend).lower()
    if name == "cpp":
        from .cpp import CppSolver
        return CppSolver(beta, info_A, info_B, info_C)
    if name == "numpy":
        from .numpy import NumpySolver
        return NumpySolver(beta, info_A, info_B, info_C)
    if name == "cupy":
        from .cupy import CupySolver
        return CupySolver(beta, info_A, info_B, info_C)
    raise ValueError(
        f"Unknown backend {backend!r}; expected one of 'cpp', 'numpy', 'cupy'."
    )
```

Now wire the TOML parser. In `python/package/chiq/bse_toml.py`, in the function that
builds the `[chiq_main]` tool dict (search for where `work_dir` / `num_wf` are read),
add a `backend` key:

```python
# in the chiq_main section parsing of bse_toml.load_params_from_toml
dict_tool["backend"] = dict_main.get("backend", "cpp")
if str(dict_tool["backend"]).lower() not in ("cpp", "numpy", "cupy"):
    raise ValueError(
        f"[chiq_main] backend must be 'cpp', 'numpy', or 'cupy', got "
        f"{dict_tool['backend']!r}"
    )
```

(Adapt the local variable name to the existing code: `dict_main` is whatever dict holds
the `[chiq_main]` table; `dict_tool` is the returned tool dict. If `num_wf` is read as
`dict_tool["num_wf"] = dict_main.get("num_wf", None)`, mirror that exact pattern.)

Now route `chiq_main.py`. Replace the direct construction at `chiq_main.py:551`:

```python
# OLD:
# solver = bse_solver.BSESolver(chi_worker.beta, matinfo_A, matinfo_B, matinfo_C)
# NEW:
from chiq.solver import get_solver  # add near the top imports; remove `import bse_solver`
...
solver = get_solver(dict_tool["backend"], chi_worker.beta, matinfo_A, matinfo_B, matinfo_C)
```

The subsequent `solver.set(...)`, `solver.calc(calc_type)`, `solver.get(...)` calls are
unchanged (the facade mirrors the pybind API). Read `backend` from the already-parsed
`dict_tool` in `main()` (it is available where `n_iw = dict_tool["num_wf"]` is read):

```python
backend = dict_tool["backend"]
```

and thread it into the loop (or read `dict_tool["backend"]` at the `get_solver` call site).

- [ ] **Step 4: Run tests to verify they pass**

Run: `python -m pytest tests/python/non-mpi/solver/test_factory_cpp.py -v`
Expected: PASS (2 tests).
Run the full existing suite to confirm the refactor is behavior-preserving:
Run: `python -m pytest tests/python/non-mpi -v`
Expected: PASS (all previously-passing tests still pass; default backend cpp).

- [ ] **Step 5: Commit**

```bash
git add python/package/chiq/solver/base.py python/package/chiq/solver/cpp.py python/package/chiq/solver/__init__.py python/package/chiq/bse_toml.py python/scripts/chiq_main.py tests/python/non-mpi/solver/test_factory_cpp.py
git commit -m "feat(solver): SolverBase + get_solver factory + CppSolver; route chiq_main via backend"
```

---

### Task 8: kernels — chi0, and the NumpySolver facade

**Files:**
- Create: `python/package/chiq/solver/kernels.py`
- Create: `python/package/chiq/solver/numpy.py`
- Test: `tests/python/non-mpi/solver/test_numpy_chi0.py`, `tests/python/non-mpi/solver/test_numpy_facade.py`

**Interfaces:**
- Consumes: `layout`, `base.SET_LAYOUT`/`CALC_TYPES`.
- Produces: `kernels.calc_chi0(state, beta, nb, nw, n_in) -> dict[name -> bm]`; `NumpySolver(SolverBase)` with `set/calc/get`, storing inputs in `self._in` (normalized complex128) and outputs in `self._out`. `set` unknown name → `RuntimeError("Invalid type '<name>'")`; missing input on `calc` → `RuntimeError("'<internal>' must be set before calling calc.")` with `X0_loc`→`'X0_Loc'`; `get` of un-computed valid name → `{}`; `get` unknown → `RuntimeError`.

- [ ] **Step 1: Write the failing test**

```python
# tests/python/non-mpi/solver/test_numpy_facade.py
import numpy as np
import pytest
from chiq.solver import get_solver

def _info(nb, nw, n_in):
    return ([n_in * nw] * nb, [n_in] * (nb * nw), [n_in] * nb)

def test_set_unknown_name_raises():
    s = get_solver("numpy", 10.0, *_info(1, 2, 1))
    with pytest.raises(RuntimeError, match="Invalid type"):
        s.set({(0, 0): np.zeros((1, 1), complex)}, "nope")

def test_get_before_calc_returns_empty_dict():
    s = get_solver("numpy", 10.0, *_info(1, 2, 1))
    assert s.get("chi0_q") == {}

def test_missing_input_reports_capital_X0_Loc():
    s = get_solver("numpy", 10.0, *_info(1, 2, 1))
    s.set({(0, 0): np.zeros((1, 1), complex)}, "chi0_loc")
    s.set({(0, 0): np.zeros((1, 1), complex), (1, 1): np.zeros((1, 1), complex)}, "X0_q")
    with pytest.raises(RuntimeError, match="'X0_Loc' must be set"):
        s.calc("chi0")
```

```python
# tests/python/non-mpi/solver/test_numpy_chi0.py
import numpy as np
from chiq.solver import get_solver

def test_chi0_formula_one_block():
    # nb=1, nw=1, n_in=1, beta given; chi0_q = chi0_loc + (1/beta)*(X0_q - X0_loc)
    beta = 4.0
    iA, iB, iC = [1], [1], [1]
    s = get_solver("numpy", beta, iA, iB, iC)
    s.set({(0, 0): np.array([[7.0]], complex)}, "chi0_loc")
    s.set({(0, 0): np.array([[2.0]], complex)}, "X0_loc")
    s.set({(0, 0): np.array([[10.0]], complex)}, "X0_q")
    s.calc("chi0")
    out = s.get("chi0_q")
    assert np.allclose(out[(0, 0)], [[7.0 + (10.0 - 2.0) / beta]])
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `python -m pytest tests/python/non-mpi/solver/test_numpy_facade.py tests/python/non-mpi/solver/test_numpy_chi0.py -v`
Expected: FAIL (`chiq.solver.numpy` / `kernels` undefined).

- [ ] **Step 3: Write minimal implementation**

```python
# python/package/chiq/solver/kernels.py
"""The five susceptibility calculations on block matrices (layout ops).

Each returns a dict of named output block matrices. Operation order mirrors
src/bse.hpp verbatim (explicit block_inverse then matmul) so the results
agree with the C++ backend within tolerance.
"""

from . import layout as L


def calc_chi0(state, beta, nb, nw, n_in):
    info_B = [n_in] * (nb * nw)  # noqa: F841 (documents B dims)
    diff = L.sub(state["X0_q"], state["X0_loc"])
    chi0_q = L.add(state["chi0_loc"], L.scale(1.0 / beta, L.sumfreq_b(diff, nb, nw, n_in)))
    return {"chi0_q": chi0_q}
```

```python
# python/package/chiq/solver/numpy.py
"""NumPy backend: runs the kernels on complex128 block matrices."""

import numpy as np

from . import kernels
from .base import SolverBase, SET_LAYOUT, CALC_TYPES, GET_NAMES
from .layout import parse_matrix_info

# internal require-names used in the C++ error strings (bse.hpp)
_REQUIRE_NAME = {"X0_loc": "X0_Loc"}  # others equal their set name

# calc_type -> list of required set names (from bse.hpp require() calls)
_REQUIRED = {
    "chi0": ["X0_q", "X0_loc", "chi0_loc"],
    "bse": ["X0_q", "X0_loc", "X_loc", "chi_loc"],
    "rpa": ["gamma0", "chi0_q"],
    "rrpa": ["chi0_loc", "chi_loc", "chi0_q"],
    "scl": ["X0_q", "X0_loc", "Phi", "Phi_sum"],
}


class NumpySolver(SolverBase):
    def __init__(self, beta, info_A, info_B, info_C):
        super().__init__(beta, info_A, info_B, info_C)
        self.nb, self.nw, self.n_in = parse_matrix_info(info_A, info_B, info_C)
        self._in = {}
        self._out = {}

    def set(self, bm, name):
        if name not in SET_LAYOUT:
            raise RuntimeError(f"Invalid type '{name}'")
        self._in[name] = {k: np.asarray(v, dtype=complex) for k, v in bm.items()}

    def _require(self, calc_type):
        for name in _REQUIRED[calc_type]:
            if name not in self._in:
                internal = _REQUIRE_NAME.get(name, name)
                raise RuntimeError(f"'{internal}' must be set before calling calc.")

    def calc(self, calc_type):
        ct = str(calc_type).lower()
        if ct not in CALC_TYPES:
            raise RuntimeError(f"Invalid type '{calc_type}'")
        self._require(ct)
        fn = {
            "chi0": kernels.calc_chi0,
            # bse/rpa/rrpa/scl added in later tasks
        }[ct]
        result = fn(self._in, self.beta, self.nb, self.nw, self.n_in)
        self._out.update(result)
        # I_q and I_q_scl are one slot in C++ (both GetIq()): whichever
        # irreducible vertex was computed most recently answers to both names.
        if "I_q" in result:
            self._out["I_q_scl"] = result["I_q"]
        if "I_q_scl" in result:
            self._out["I_q"] = result["I_q_scl"]

    def get(self, name):
        # unknown output name -> RuntimeError (match C++ getMatrix); valid but
        # not-yet-computed -> {}
        if name not in GET_NAMES:
            raise RuntimeError(f"Invalid type '{name}'")
        return self._out.get(name, {})
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `python -m pytest tests/python/non-mpi/solver/test_numpy_facade.py tests/python/non-mpi/solver/test_numpy_chi0.py -v`
Expected: PASS (4 tests).

- [ ] **Step 5: Commit**

```bash
git add python/package/chiq/solver/kernels.py python/package/chiq/solver/numpy.py tests/python/non-mpi/solver/test_numpy_facade.py tests/python/non-mpi/solver/test_numpy_chi0.py
git commit -m "feat(solver): NumpySolver facade + chi0 kernel"
```

---

### Task 9: kernels — bse

**Files:**
- Modify: `python/package/chiq/solver/kernels.py`
- Modify: `python/package/chiq/solver/numpy.py` (register `bse`)
- Test: `tests/python/non-mpi/solver/test_numpy_bse.py`

**Interfaces:**
- Consumes: `layout`, `calc_chi0`.
- Produces: `kernels.calc_bse(state, beta, nb, nw, n_in) -> {"chi_q":..., "I_q":...}`.

- [ ] **Step 1: Write the failing test**

```python
# tests/python/non-mpi/solver/test_numpy_bse.py
import numpy as np
from chiq.solver import get_solver
from chiq.solver.kernels import calc_bse

def _info(nb, nw, n_in):
    return ([n_in * nw] * nb, [n_in] * (nb * nw), [n_in] * nb)

def test_iq_and_iq_scl_share_one_slot():
    # After a bse calc, get("I_q_scl") returns the same value as get("I_q").
    nb, nw, n_in = 1, 2, 1
    s = get_solver("numpy", 5.0, *_info(nb, nw, n_in))
    X0 = {(0, 0): np.array([[1.0]], complex), (1, 1): np.array([[2.0]], complex)}
    s.set(X0, "X0_q")
    s.set(X0, "X0_loc")
    s.set({(0, 0): np.eye(2, dtype=complex)}, "X_loc")
    s.set({(0, 0): np.array([[3.0]], complex)}, "chi_loc")
    s.calc("bse")
    assert s.get("I_q_scl") == s.get("I_q")

def test_bse_reduces_to_chi_loc_when_X0q_equals_X0loc():
    # If X0_q == X0_loc then P_q = 0, Z_q = 0, so chi_q == chi_loc.
    nb, nw, n_in, beta = 1, 2, 1, 5.0
    X0 = {(0, 0): np.array([[1.0]], complex), (1, 1): np.array([[2.0]], complex)}
    X_loc = {(0, 0): np.array([[1.0, 0.0], [0.0, 1.0]], complex)}  # A-type 2x2
    chi_loc = {(0, 0): np.array([[3.0]], complex)}
    state = {"X0_q": X0, "X0_loc": X0, "X_loc": X_loc, "chi_loc": chi_loc}
    out = calc_bse(state, beta, nb, nw, n_in)
    assert np.allclose(out["chi_q"][(0, 0)], [[3.0]])
    # I_q = chi_loc^-1 - chi_q^-1 = 0 here
    assert np.allclose(out["I_q"][(0, 0)], [[0.0]])
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python -m pytest tests/python/non-mpi/solver/test_numpy_bse.py -v`
Expected: FAIL (`calc_bse` undefined).

- [ ] **Step 3: Write minimal implementation**

```python
# append to python/package/chiq/solver/kernels.py

def calc_bse(state, beta, nb, nw, n_in):
    info_A = [n_in * nw] * nb
    info_B = [n_in] * (nb * nw)  # noqa: F841
    info_C = [n_in] * nb  # noqa: F841
    X0_loc = state["X0_loc"]
    X0_q = state["X0_q"]
    X_loc = state["X_loc"]
    chi_loc = state["chi_loc"]

    Pq = L.sub(L.block_inverse(X0_loc, [n_in] * (nb * nw)),
               L.block_inverse(X0_q, [n_in] * (nb * nw)))            # B-type
    Xloc_Pq = L.prod_ab(X_loc, Pq, nb, nw, n_in)                     # A-type
    ident_A = L.identity(info_A)
    inner = L.sub(L.block_inverse(L.sub(ident_A, Xloc_Pq), info_A), ident_A)  # A-type
    Zq = L.matmul(inner, X_loc, info_A)                             # A-type
    chi_q = L.add(chi_loc, L.scale(1.0 / beta, L.sumfreq_a(Zq, nb, nw, n_in)))  # C-type
    Iq = L.sub(L.block_inverse(chi_loc, [n_in] * nb),
               L.block_inverse(chi_q, [n_in] * nb))                  # C-type
    return {"chi_q": chi_q, "I_q": Iq}
```

Register it in `numpy.py`'s dispatch dict:

```python
# in NumpySolver.calc, extend the dispatch:
fn = {
    "chi0": kernels.calc_chi0,
    "bse": kernels.calc_bse,
}[ct]
```

- [ ] **Step 4: Run test to verify it passes**

Run: `python -m pytest tests/python/non-mpi/solver/test_numpy_bse.py -v`
Expected: PASS (1 test).

- [ ] **Step 5: Commit**

```bash
git add python/package/chiq/solver/kernels.py python/package/chiq/solver/numpy.py tests/python/non-mpi/solver/test_numpy_bse.py
git commit -m "feat(solver): bse kernel"
```

---

### Task 10: kernels — rpa and rrpa

**Files:**
- Modify: `python/package/chiq/solver/kernels.py`
- Modify: `python/package/chiq/solver/numpy.py`
- Test: `tests/python/non-mpi/solver/test_numpy_rpa.py`

**Interfaces:**
- Consumes: `layout`.
- Produces: `kernels.calc_rpa(state, beta, nb, nw, n_in) -> {"chi_q_rpa":...}`; `kernels.calc_rrpa(...) -> {"chi_q_rrpa":...}`.

- [ ] **Step 1: Write the failing test**

```python
# tests/python/non-mpi/solver/test_numpy_rpa.py
import numpy as np
from chiq.solver.kernels import calc_rpa, calc_rrpa

def test_rpa_zero_gamma_gives_chi0():
    nb, nw, n_in, beta = 1, 1, 1, 1.0
    chi0_q = {(0, 0): np.array([[2.0]], complex)}
    gamma0 = {(0, 0): np.array([[0.0]], complex)}
    out = calc_rpa({"chi0_q": chi0_q, "gamma0": gamma0}, beta, nb, nw, n_in)
    assert np.allclose(out["chi_q_rpa"][(0, 0)], [[2.0]])

def test_rrpa_when_chi0loc_equals_chiloc_gives_chi0q():
    # U_eff = chi0_loc^-1 - chi_loc^-1 = 0 -> chi_q_rrpa = chi0_q
    nb, nw, n_in, beta = 1, 1, 1, 1.0
    same = {(0, 0): np.array([[3.0]], complex)}
    chi0_q = {(0, 0): np.array([[2.0]], complex)}
    out = calc_rrpa({"chi0_loc": same, "chi_loc": same, "chi0_q": chi0_q}, beta, nb, nw, n_in)
    assert np.allclose(out["chi_q_rrpa"][(0, 0)], [[2.0]])
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python -m pytest tests/python/non-mpi/solver/test_numpy_rpa.py -v`
Expected: FAIL (`calc_rpa`/`calc_rrpa` undefined).

- [ ] **Step 3: Write minimal implementation**

```python
# append to python/package/chiq/solver/kernels.py

def calc_rpa(state, beta, nb, nw, n_in):
    info_C = [n_in] * nb
    chi0_q = state["chi0_q"]
    gamma0 = state["gamma0"]
    m = L.sub(L.identity(info_C), L.matmul(chi0_q, gamma0, info_C))
    chi_q_rpa = L.matmul(L.block_inverse(m, info_C), chi0_q, info_C)
    return {"chi_q_rpa": chi_q_rpa}


def calc_rrpa(state, beta, nb, nw, n_in):
    info_C = [n_in] * nb
    chi0_loc = state["chi0_loc"]
    chi_loc = state["chi_loc"]
    chi0_q = state["chi0_q"]
    Ueff = L.sub(L.block_inverse(chi0_loc, info_C), L.block_inverse(chi_loc, info_C))
    m = L.sub(L.identity(info_C), L.matmul(chi0_q, Ueff, info_C))
    chi_q_rrpa = L.matmul(L.block_inverse(m, info_C), chi0_q, info_C)
    return {"chi_q_rrpa": chi_q_rrpa}
```

Register both in `numpy.py`:

```python
fn = {
    "chi0": kernels.calc_chi0,
    "bse": kernels.calc_bse,
    "rpa": kernels.calc_rpa,
    "rrpa": kernels.calc_rrpa,
}[ct]
```

- [ ] **Step 4: Run test to verify it passes**

Run: `python -m pytest tests/python/non-mpi/solver/test_numpy_rpa.py -v`
Expected: PASS (2 tests).

- [ ] **Step 5: Commit**

```bash
git add python/package/chiq/solver/kernels.py python/package/chiq/solver/numpy.py tests/python/non-mpi/solver/test_numpy_rpa.py
git commit -m "feat(solver): rpa and rrpa kernels"
```

---

### Task 11: kernels — scl

**Files:**
- Modify: `python/package/chiq/solver/kernels.py`
- Modify: `python/package/chiq/solver/numpy.py`
- Test: `tests/python/non-mpi/solver/test_numpy_scl.py`

**Interfaces:**
- Consumes: `layout`.
- Produces: `kernels.calc_scl(state, beta, nb, nw, n_in) -> {"chi_q_scl":..., "I_q_scl":...}`. Note: `Phi` (B-type) and `Phi_sum` (C-type) are provided by the driver via `set`; the kernel does not compute the matrix sqrt.

- [ ] **Step 1: Write the failing test**

```python
# tests/python/non-mpi/solver/test_numpy_scl.py
import numpy as np
from chiq.solver.kernels import calc_scl

def test_scl_when_X0q_equals_X0loc_gives_phisum_squared():
    # Lambda_q = X0_loc^-1 - X0_q^-1 = 0 -> lambda_q = 0
    # chi_q_scl = Phi_sum * I^-1 * Phi_sum = Phi_sum @ Phi_sum
    nb, nw, n_in, beta = 1, 2, 1, 1.0
    X0 = {(0, 0): np.array([[1.0]], complex), (1, 1): np.array([[2.0]], complex)}
    Phi = {(0, 0): np.array([[1.0]], complex), (1, 1): np.array([[1.0]], complex)}  # B-type
    Phi_sum = {(0, 0): np.array([[3.0]], complex)}  # C-type
    out = calc_scl({"X0_q": X0, "X0_loc": X0, "Phi": Phi, "Phi_sum": Phi_sum}, beta, nb, nw, n_in)
    assert np.allclose(out["chi_q_scl"][(0, 0)], [[9.0]])   # 3*3
    assert np.allclose(out["I_q_scl"][(0, 0)], [[0.0]])     # Phi_sum^-1 * 0 * Phi_sum^-1
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python -m pytest tests/python/non-mpi/solver/test_numpy_scl.py -v`
Expected: FAIL (`calc_scl` undefined).

- [ ] **Step 3: Write minimal implementation**

```python
# append to python/package/chiq/solver/kernels.py

def calc_scl(state, beta, nb, nw, n_in):
    info_B = [n_in] * (nb * nw)
    info_C = [n_in] * nb
    X0_loc = state["X0_loc"]
    X0_q = state["X0_q"]
    Phi = state["Phi"]
    Phi_sum = state["Phi_sum"]

    Lambdaq = L.sub(L.block_inverse(X0_loc, info_B), L.block_inverse(X0_q, info_B))  # B
    lambda_q = L.sumfreq_b(L.matmul(L.matmul(Phi, Lambdaq, info_B), Phi, info_B),
                           nb, nw, n_in)                                            # C
    m = L.sub(L.identity(info_C), lambda_q)
    chi_q_scl = L.matmul(L.matmul(Phi_sum, L.block_inverse(m, info_C), info_C),
                         Phi_sum, info_C)
    phisum_inv = L.block_inverse(Phi_sum, info_C)
    Iq = L.matmul(L.matmul(phisum_inv, lambda_q, info_C), phisum_inv, info_C)
    return {"chi_q_scl": chi_q_scl, "I_q_scl": Iq}
```

Register in `numpy.py`:

```python
fn = {
    "chi0": kernels.calc_chi0,
    "bse": kernels.calc_bse,
    "rpa": kernels.calc_rpa,
    "rrpa": kernels.calc_rrpa,
    "scl": kernels.calc_scl,
}[ct]
```

- [ ] **Step 4: Run test to verify it passes**

Run: `python -m pytest tests/python/non-mpi/solver/test_numpy_scl.py -v`
Expected: PASS (1 test).

- [ ] **Step 5: Commit**

```bash
git add python/package/chiq/solver/kernels.py python/package/chiq/solver/numpy.py tests/python/non-mpi/solver/test_numpy_scl.py
git commit -m "feat(solver): scl kernel"
```

---

### Task 12: backend-agreement tests over the integration fixtures

**Files:**
- Create: `tests/python/non-mpi/solver/test_backend_agreement.py`
- Test: itself.

**Interfaces:**
- Consumes: `get_solver`, existing fixture HDF5 inputs and the `chi_*` worker classes in `chiq_main` (imported for building the per-(ω,q) block-matrix dicts), OR — simpler and preferred — construct random well-conditioned inputs of the shapes each calc needs and compare `cpp` vs `numpy` directly. Use the direct approach to avoid coupling to the driver's HDF5 layout.

- [ ] **Step 1: Write the failing test**

```python
# tests/python/non-mpi/solver/test_backend_agreement.py
import numpy as np
import pytest
from chiq.solver import get_solver

RNG = np.random.default_rng(0)

def _info(nb, nw, n_in):
    return ([n_in * nw] * nb, [n_in] * (nb * nw), [n_in] * nb)

def _rand_c(shape):
    return (RNG.standard_normal(shape) + 1j * RNG.standard_normal(shape))

def _well_conditioned(n):
    # diagonally dominant -> invertible, bounded condition number
    m = _rand_c((n, n))
    return m + n * np.eye(n)

def _assert_agree(a, b):
    for k in set(a) | set(b):
        if k in a and k in b:
            va, vb = a[k], b[k]
        elif k in a:
            va, vb = a[k], np.zeros_like(a[k])
        else:
            va, vb = np.zeros_like(b[k]), b[k]
        # NaN / signed-inf position masks, then elementwise mixed tolerance.
        # np.isposinf/np.isneginf reject complex dtype, so check real/imag parts
        # separately (a test-harness detail, not a loosening of the metric).
        va_r, va_i = np.real(va), np.imag(va)
        vb_r, vb_i = np.real(vb), np.imag(vb)
        assert np.array_equal(np.isnan(va_r), np.isnan(vb_r))
        assert np.array_equal(np.isnan(va_i), np.isnan(vb_i))
        assert np.array_equal(np.isposinf(va_r), np.isposinf(vb_r))
        assert np.array_equal(np.isposinf(va_i), np.isposinf(vb_i))
        assert np.array_equal(np.isneginf(va_r), np.isneginf(vb_r))
        assert np.array_equal(np.isneginf(va_i), np.isneginf(vb_i))
        np.testing.assert_allclose(va, vb, rtol=1e-8, atol=1e-10, equal_nan=True)

def _c_blocks(nb, n_in):
    return {(i, j): _rand_c((n_in, n_in)) for i in range(nb) for j in range(nb)}

def _c_blocks_wc(nb, n_in):
    # C-type with well-conditioned diagonal blocks (invertible)
    d = {(i, i): _well_conditioned(n_in) for i in range(nb)}
    return d

def _b_blocks_diag(nb, nw, n_in):
    return {(iw * nb + b, iw * nb + b2): _rand_c((n_in, n_in))
            for iw in range(nw) for b in range(nb) for b2 in range(nb)}

def _b_blocks_diag_wc(nb, nw, n_in):
    # B-type well-conditioned in each same-iw super-block (diagonal blocks dominant)
    out = {}
    for iw in range(nw):
        for b in range(nb):
            for b2 in range(nb):
                m = _rand_c((n_in, n_in))
                if b == b2:
                    m = m + (nb * n_in) * np.eye(n_in)
                out[(iw * nb + b, iw * nb + b2)] = m
    return out

def _a_blocks(nb, nw, n_in):
    return {(i, j): _rand_c((n_in * nw, n_in * nw)) for i in range(nb) for j in range(nb)}

@pytest.mark.parametrize("nb,nw,n_in", [(1, 2, 1), (2, 2, 2), (2, 3, 2)])
def test_chi0_agreement(nb, nw, n_in):
    iA, iB, iC = _info(nb, nw, n_in)
    inputs = {
        "chi0_loc": {(i, j): _rand_c((n_in, n_in)) for i in range(nb) for j in range(nb)},
        "X0_loc": {(iw * nb + b, iw * nb + b2): _rand_c((n_in, n_in))
                   for iw in range(nw) for b in range(nb) for b2 in range(nb)},
    }
    inputs["X0_q"] = {k: _rand_c((n_in, n_in)) for k in inputs["X0_loc"]}
    outs = {}
    for backend in ("cpp", "numpy"):
        s = get_solver(backend, 12.0, iA, iB, iC)
        for name, bm in inputs.items():
            s.set(bm, name)
        s.calc("chi0")
        outs[backend] = s.get("chi0_q")
    _assert_agree(outs["cpp"], outs["numpy"])

@pytest.mark.parametrize("nb,nw,n_in", [(1, 2, 1), (2, 2, 2)])
def test_rpa_agreement(nb, nw, n_in):
    iA, iB, iC = _info(nb, nw, n_in)
    chi0_q = {(i, j): _rand_c((n_in, n_in)) for i in range(nb) for j in range(nb)}
    gamma0 = {(i, i): 0.01 * _rand_c((n_in, n_in)) for i in range(nb)}
    outs = {}
    for backend in ("cpp", "numpy"):
        s = get_solver(backend, 12.0, iA, iB, iC)
        s.set(chi0_q, "chi0_q")
        s.set(gamma0, "gamma0")
        s.calc("rpa")
        outs[backend] = s.get("chi_q_rpa")
    _assert_agree(outs["cpp"], outs["numpy"])

@pytest.mark.parametrize("nb,nw,n_in", [(1, 2, 1), (2, 2, 2)])
def test_rrpa_agreement(nb, nw, n_in):
    iA, iB, iC = _info(nb, nw, n_in)
    inputs = {
        "chi0_loc": _c_blocks_wc(nb, n_in),
        "chi_loc": _c_blocks_wc(nb, n_in),
        "chi0_q": _c_blocks(nb, n_in),
    }
    outs = {}
    for backend in ("cpp", "numpy"):
        s = get_solver(backend, 12.0, iA, iB, iC)
        for name, bm in inputs.items():
            s.set(bm, name)
        s.calc("rrpa")
        outs[backend] = s.get("chi_q_rrpa")
    _assert_agree(outs["cpp"], outs["numpy"])

@pytest.mark.parametrize("nb,nw,n_in", [(1, 2, 1), (2, 2, 2)])
def test_bse_agreement(nb, nw, n_in):
    iA, iB, iC = _info(nb, nw, n_in)
    X0_loc = _b_blocks_diag_wc(nb, nw, n_in)
    inputs = {
        "X0_loc": X0_loc,
        "X0_q": _b_blocks_diag_wc(nb, nw, n_in),
        "X_loc": _a_blocks(nb, nw, n_in),
        "chi_loc": _c_blocks_wc(nb, n_in),
    }
    outs = {}
    for backend in ("cpp", "numpy"):
        s = get_solver(backend, 12.0, iA, iB, iC)
        for name, bm in inputs.items():
            s.set(bm, name)
        s.calc("bse")
        outs[backend] = s.get("chi_q")
    _assert_agree(outs["cpp"], outs["numpy"])

@pytest.mark.parametrize("nb,nw,n_in", [(1, 2, 1), (2, 2, 2)])
def test_scl_agreement(nb, nw, n_in):
    iA, iB, iC = _info(nb, nw, n_in)
    inputs = {
        "X0_loc": _b_blocks_diag_wc(nb, nw, n_in),
        "X0_q": _b_blocks_diag_wc(nb, nw, n_in),
        "Phi": _b_blocks_diag(nb, nw, n_in),
        "Phi_sum": _c_blocks_wc(nb, n_in),
    }
    outs = {}
    for backend in ("cpp", "numpy"):
        s = get_solver(backend, 12.0, iA, iB, iC)
        for name, bm in inputs.items():
            s.set(bm, name)
        s.calc("scl")
        outs[backend] = s.get("chi_q_scl")
    _assert_agree(outs["cpp"], outs["numpy"])
```

- [ ] **Step 2: Run test to verify it fails or errors**

Run: `python -m pytest tests/python/non-mpi/solver/test_backend_agreement.py -v` (needs built `bse_solver`).
Expected: initially FAIL if any kernel disagrees with C++ — investigate and fix the layout/kernel until agreement holds. (Common causes: a transpose in `sumfreq_a` reshape, or a wrong operand order in a `matmul`.)

- [ ] **Step 3: Fix any disagreements**

If a test fails, compare a single failing block between backends, localize to one layout
op with a targeted unit test, fix `layout.py`/`kernels.py`, and re-run. Do not loosen the
tolerance; the metric is fixed by the Global Constraints.

- [ ] **Step 4: Run test to verify it passes**

Run: `python -m pytest tests/python/non-mpi/solver/test_backend_agreement.py -v`
Expected: PASS (all parametrizations).

- [ ] **Step 5: Commit**

```bash
git add tests/python/non-mpi/solver/test_backend_agreement.py
git commit -m "test(solver): cpp vs numpy backend agreement across calc types"
```

---

### Task 13: CupySolver skeleton + factory wiring

**Files:**
- Create: `python/package/chiq/solver/cupy.py`
- Test: `tests/python/non-mpi/solver/test_cupy_skeleton.py`

**Interfaces:**
- Consumes: `SolverBase`.
- Produces: `CupySolver(SolverBase)` whose constructor raises a clear error (cupy not enabled / not installed). No numerical execution this phase.

- [ ] **Step 1: Write the failing test**

```python
# tests/python/non-mpi/solver/test_cupy_skeleton.py
import pytest
from chiq.solver import get_solver

def test_cupy_backend_not_enabled():
    with pytest.raises((NotImplementedError, ImportError, RuntimeError)) as exc:
        get_solver("cupy", 10.0, [2], [1, 1], [1])
    assert "cupy" in str(exc.value).lower()
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python -m pytest tests/python/non-mpi/solver/test_cupy_skeleton.py -v`
Expected: FAIL (`chiq.solver.cupy` undefined -> factory ImportError with wrong message, or ModuleNotFoundError).

- [ ] **Step 3: Write minimal implementation**

```python
# python/package/chiq/solver/cupy.py
"""CuPy/GPU backend skeleton. Enabled in a later GPU phase."""

from .base import SolverBase


class CupySolver(SolverBase):
    def __init__(self, beta, info_A, info_B, info_C):
        super().__init__(beta, info_A, info_B, info_C)
        raise NotImplementedError(
            "cupy backend not yet enabled (planned for the GPU phase)."
        )
```

- [ ] **Step 4: Run test to verify it passes**

Run: `python -m pytest tests/python/non-mpi/solver/test_cupy_skeleton.py -v`
Expected: PASS (1 test).

- [ ] **Step 5: Commit**

```bash
git add python/package/chiq/solver/cupy.py tests/python/non-mpi/solver/test_cupy_skeleton.py
git commit -m "feat(solver): CupySolver skeleton raising a clear not-enabled error"
```

---

## Final verification

- [ ] Run the whole non-MPI suite (with the built `bse_solver` importable):
  `python -m pytest tests/python/non-mpi -v` — all green, including the pre-existing fixtures (proving the Task 7 refactor is behavior-preserving) and the new `solver/` tests.
- [ ] Confirm default behavior unchanged: run an example end-to-end (`examples/square_bse/run.sh` with `NPROCS=1`) and diff outputs against a pre-change run — identical.
- [ ] Confirm the `numpy` backend runs end-to-end: set `[chiq_main] backend = "numpy"` in a fixture `bse.toml` and run `chiq_main.py` on it; results match the `cpp` run within tolerance.

## Out of scope (separate plan)

pip packaging (§5 of the spec): `pyproject.toml` + scikit-build-core, renaming the pybind
module/target to `_bse_solver`, `console_scripts` entry points, `bse`/`bse_solver`
compatibility shims, wheel/sdist install tests, Python-floor bump to 3.10, and the
`toml`/`scipy` dependency declarations. That plan swaps `cpp.py`'s `import bse_solver`
for `from chiq import _bse_solver` behind the top-level shim.
