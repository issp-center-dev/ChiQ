# Design: NumPy solver backend + backend-switch + pip packaging (Phase 0)

**Date:** 2026-07-12
**Status:** Draft for review (round 3 — incorporates Codex design review rounds 1–3)
**Repo:** ChiQ (Bethe-Salpeter / DMFT momentum-dependent susceptibility solver)

## 1. Context and goal

ChiQ's heavy numerics live in a C++ BSE solver (`src/bse.hpp`, `src/block_matrix.hpp`)
exposed to Python via pybind11 (`src/bse_solver_pybind.cpp`, module `bse_solver`).
Python scripts (`python/scripts/chiq_main.py`, `chiq_post.py`) drive it: for each
(bosonic frequency ω, momentum q) they build a block matrix as a Python
`dict{(block1, block2): ndarray}`, call `solver.set(...)`, `solver.calc(type)`,
then `solver.get(...)`.

The long-term goal (agreed with the user) is **GPU acceleration** and **IR-basis
(sparse-ir) adoption**. The chosen strategy is *not* to CUDA-ify the C++, but to
introduce a **NumPy/CuPy-switchable solver backend in Python**, because the data
already lives as dict-of-ndarrays and the operations are dense block linear algebra.

**This spec (Phase 0)** is the foundation everything else builds on. Three deliverables:

1. A **backend abstraction** (`backend = "cpp" | "numpy" | "cupy"`) with the C++
   and NumPy backends implemented (CuPy is skeleton only, activated in a later phase).
2. **Numerical-agreement tests** proving the NumPy backend matches C++.
3. **pip packaging** via scikit-build-core so `pip install .` builds the C++
   extension and installs the CLI, replacing the CMake-install + `chiqvars.sh`
   PYTHONPATH workflow.

Out of scope: sparse-ir integration, actual GPU/CuPy execution and tuning, and any
change to the BSE numerical algorithm. Those are later phases that depend on this one.

## 2. Non-goals / constraints

- **Zero behavior change by default.** Default backend stays `cpp`; existing runs,
  outputs, and the `import bse` compatibility path are preserved.
- **C++ is retained**, not retired. It is a co-equal, permanent backend and the
  numerical oracle for the NumPy/CuPy backends.
- **Compatibility target = mathematical equivalence within tolerance**, not bitwise
  reproduction of C++ rounding. (This distinction is deliberate so future numerical
  improvements are not blocked by the oracle — see §4.)
- **Python floor rises to 3.10** (3.8/3.9 EOL). CI matrix becomes 3.10 … 3.13; the
  dependency/version matrix is pinned in §5.4.
- No unrelated refactors; targeted fixes only where they serve this work.

## 3. Architecture: backend abstraction

### 3.1 Component decomposition

Three separated concerns (Codex alternative "separate layout / state / kernels"):

- **`chiq.solver.layout`** — the *only* owner of block-dict ↔ dense mapping and the
  A/B/C conventions (§3.4). Pure functions, no solver state. `chiq_main._reformat_data`
  and the `_convert_block_matrix` reshapes move here; `chiq_main.py` keeps no
  numerical/indexing logic.
- **`chiq.solver` kernels** — the five calculations (§3.5) written against a narrow
  array/linalg protocol `xp` (an injected module: `numpy` or `cupy`) plus a tiny
  `linalg` adapter (`solve`, `inv`). One shared kernel implementation; the backend is
  chosen by which `xp`/adapter is injected — **not** by subclassing (Codex should-fix).
- **`chiq.solver` facade** — the stateful `set`/`calc`/`get` triad (§3.2, §3.3),
  a thin compatibility layer over the kernels.

### 3.2 Public interface (facade)

```
class SolverBase:
    def __init__(self, beta, matrix_info_A, matrix_info_B, matrix_info_C): ...
    def set(self, block_matrix_dict, name: str) -> None
    def calc(self, calc_type: str) -> None      # "chi0"|"bse"|"rpa"|"rrpa"|"scl"
    def get(self, name: str) -> dict             # {(b1,b2): ndarray}
```

Backends:

- `CppSolver` — thin wrapper over `bse_solver.BSESolver`. Behavior identical to today
  (delegates set/calc/get verbatim; complex128; same error strings).
- `NumpySolver` — kernels with `xp = numpy`, `linalg = numpy.linalg`.
- `CupySolver` — same kernels with `xp = cupy`; host→device on `set`, device→host on
  `get`. **Skeleton only** this phase: constructed via the factory but raises a clear
  `NotImplementedError("cupy backend not yet enabled")` / a clear error if cupy is
  absent. No CI execution (no GPU).

`get_solver(backend, beta, info_A, info_B, info_C) -> SolverBase` is the factory.
`chiq_main.py` obtains its solver only through it.

### 3.3 State machine and per-calc contract (normative)

**Set names → matrix-info type** (from `src/bse_solver_pybind.cpp:46-100`):

| name       | layout | notes                         |
|------------|--------|-------------------------------|
| `X_loc`    | A      | local full two-particle       |
| `X0_loc`   | B      | local bare                    |
| `X0_q`     | B      | q-dependent bare              |
| `Phi`      | B      | SCL vertex                    |
| `chi_loc`  | C      | local physical susceptibility |
| `chi0_loc` | C      | local bare physical           |
| `chi0_q`   | C      | q bare physical (RPA/RRPA in) |
| `gamma0`   | C      | bare vertex                   |
| `Phi_sum`  | C      | SCL summed vertex (√chi_loc)  |

**calc_type → required set names / produced get names** (from `bse.hpp` `require()`
calls and `chiq_main.py` worker classes `chi_chi0/chi_bse/chi_rpa/chi_rrpa/chi_scl`):

| calc   | requires (set)                    | produces (get)                 |
|--------|-----------------------------------|--------------------------------|
| `chi0` | X0_q, X0_loc, chi0_loc            | chi0_q                         |
| `bse`  | X0_q, X0_loc, X_loc, chi_loc      | chi0_q?†, chi_q, I_q           |
| `rpa`  | chi0_q, gamma0                    | chi_q_rpa                      |
| `rrpa` | chi0_q, chi0_loc, chi_loc         | chi_q_rrpa                     |
| `scl`  | X0_q, X0_loc, Phi, Phi_sum        | chi_q_scl, I_q_scl             |

† `chi0_q` for `bse` is produced by a preceding `chi0` calc in the driver, not by
`CalcChiq_BSE`; the driver's `output_q` lists reflect what is written to HDF5. The
solver contract is: `get(name)` returns the last computed result of that name and
raises `KeyError`-style errors if that calc has not run.

**Lifecycle / error semantics** (match observed C++ behavior):
- `set` is idempotent-overwrite: setting the same name replaces prior data.
- `calc` requires all inputs in its row to be set; a missing input raises
  `RuntimeError("'<name>' must be set before calling calc.")` (mirrors
  `bse.hpp` `require`).
- Repeated `calc` recomputes from currently-set inputs and overwrites outputs.
- **No output caching / no auto-invalidation** (matches C++ `BSESolver`, which holds
  inputs by pointer and recomputes on each `calc`): `get` returns the result of the
  most recent `calc`; changing an input via `set` does **not** retroactively update or
  invalidate a prior `get` result. The driver's contract is `set → calc → get` in
  order per (ω,q), so staleness cannot arise in practice; this behavior is documented
  as compatibility-required rather than "fixed".
- Unknown set/get/calc names raise `RuntimeError("Invalid type '<name>'")` (mirrors
  the pybind `else` branches). `CppSolver` inherits these verbatim; `NumpySolver`
  raises the same messages so tests are backend-agnostic.

**Ownership/aliasing:** `set` copies inputs (defensive; matches C++ which copies into
its own block_matrix). `get` returns freshly-assembled dicts of owned arrays (not
internal views). A round-trip mutation test asserts this for both backends.

**dtype:** inputs are normalized to `complex128` on `set` (the C++ binding does the
same via `dict_to_blockmatrix<complex<double>>`); `get` returns `complex128`.

### 3.4 Canonical block layout (normative)

**Every block matrix is a `dict{(vi, vj): ndarray}`** keyed by an *ordered pair of
integer vertices* `(vi, vj)`; a missing key means a **structural zero** (a coupling
that is mathematically forbidden, not merely numerically small — see §8 risk
disposition). Dimensions come solely from the three matrix-info vectors; no global
state is consulted. Let `nb` = number of C-vertices, `nw` = number of frequencies,
`n_in` = inner dimension.

**Uniform-dimension invariant (committed).** In ChiQ, `n_in` and `nw` are the same
for every block: the driver builds `matinfo_A = [n_iw*n_inner]*n_block`,
`matinfo_B = [n_inner]*(n_block*n_iw)`, `matinfo_C = [n_inner]*n_block`
(`chiq_main.py:524-526`), and `_save_X0_info` asserts identical inner structure
across all blocks (`sumk_dft_chi.py:579-580`). The layout module **requires** uniform
`n_in`/`nw` and raises on any non-uniform matrix-info. (The per-component solve in
§3.5 does not itself rely on uniformity, but the public contract does; heterogeneous
blocks are explicitly unsupported in Phase 0.)

| layout | vertex set (key domain)        | array shape per key      | info vector       |
|--------|--------------------------------|--------------------------|-------------------|
| C      | `v ∈ [0, nb)`                  | `n_in × n_in`            | `info_C[v]=n_in`  |
| B      | `v = iw*nb + b`, `iw∈[0,nw)`   | `n_in × n_in`            | `info_B[v]=n_in`  |
| A      | `v ∈ [0, nb)`                  | `(n_in*nw) × (n_in*nw)`  | `info_A[v]=n_in*nw` |

- B is structurally block-diagonal in `iw`: only keys `(iw*nb+b_i, iw*nb+b_j)`
  (same `iw`) may be present.
- A flattens its inner index as `in*nw + iw`.
- All blocks are **square**; the layout module asserts this and A/B/C mutual
  consistency (`nb, nw, n_in`) on construction, raising on mismatch.

Reductions/conversions are functions `dict → dict` with the exact key/shape contract
(reproduced from `bse.hpp`, each with a unit test):
- `Sumfreq_A` (A→C): output key `(b_i,b_j)`, shape `n_in×n_in`;
  `out[(i,j)][a,c] = Σ_{iw1,iw2} A[(i,j)][a*nw+iw1, c*nw+iw2]`.
- `Sumfreq_B` (B→C): output key `(b_i,b_j)`;
  `out[(i,j)] = Σ_iw B[(iw*nb+i, iw*nb+j)]` over keys present.
- `Convert_A2B` (A→B) / `Convert_B2A` (B→A): index remap per `bse.hpp:519-613`;
  keys and shapes follow the B/A domains above.

### 3.5 Kernel formulas (ordered, matching C++)

Identity matrices carry the layout of the operand. `β` scaling is applied on the
C-type sum as shown. `solve(A, B) ≡ A⁻¹ B` (left solve).

- **chi0** (`bse.hpp:223`): `chi0_q = chi0_loc + (1/β)·Sumfreq_B(X0_q − X0_loc)`.
- **bse** (`bse.hpp:244`):
  `Pq = inv(X0_loc) − inv(X0_q)` (B-type, explicit inverse — Pq is consumed as data);
  `XlocPq = ProdAB(X_loc, Pq)` (A-type);
  `Zq = solve(I − XlocPq, X_loc) − X_loc` (A-type; `([1−XlocPq]⁻¹ − 1)·X_loc`);
  `chi_q = chi_loc + (1/β)·Sumfreq_A(Zq)`;
  `I_q = inv(chi_loc) − inv(chi_q)` (explicit inverses — outputs).
- **rpa** (`bse.hpp:281`): `chi_q_rpa = solve(I − chi0_q·gamma0, chi0_q)` (C-type).
- **rrpa** (`bse.hpp:301`): `Ueff = inv(chi0_loc) − inv(chi_loc)`;
  `chi_q_rrpa = solve(I − chi0_q·Ueff, chi0_q)`.
- **scl** (`bse.hpp:326`): `Lambdaq = inv(X0_loc) − inv(X0_q)` (B);
  `lambda_q = Sumfreq_B(Phi·Lambdaq·Phi)` (C);
  `chi_q_scl = Phi_sum · solve(I − lambda_q, Phi_sum)`;
  `I_q_scl = inv(Phi_sum)·lambda_q·inv(Phi_sum)`.

**Block-structure-preserving inverse/solve + fill-in rule (resolves dense-vs-groupwise
and fill-in concerns).** `inv`/`solve` are performed **per connected component** of
the block-sparsity graph, reproducing `block_matrix.hpp:541` (`gather_blocks`):
assemble each component into a dense super-block, invert/solve it, scatter back.
Component derivation and fill-in match C++ exactly:
- **Components** are connected components under the undirected edge
  `connected(i,j) = exists(i,j) or exists(j,i)` (`block_matrix.hpp:547`). For binary
  `solve(A, B)`, components are taken from `A` (the matrix being inverted); `B` must
  share `A`'s block structure (C++ enforces this in `operator*` /
  `has_same_block_structure`).
- **Fill-in:** the inverted super-block is dense, and C++ **emits every vertex pair
  `(vi, vj)` within the component** (`block_matrix.hpp:590-602` assigns all `bl,bl2`),
  including pairs absent from the input. The NumPy `inv`/`solve` reproduce this: output
  keys = all intra-component pairs; **numerical zeros are retained** (not pruned).
  Cross-component pairs remain absent (structural zero).
- Products (`*`) emit keys **symbolically**, not by numerical value: C++ `operator*`
  (`block_matrix.hpp:639-650`) sets a key `(i,j)` iff **∃ k with both `(i,k)` and
  `(k,j)` present** — the decision is `exists(i,k) && exists(k,j)`, independent of the
  numerical result. The NumPy `*` reproduces this exact key-path rule, so the output
  sparsity graph is backend-identical regardless of floating-point cancellation.
  Differences (`-`) are entrywise over the union of present keys. Downstream
  `Sumfreq_*`/`Convert_*`/`get` consume whatever keys are present. No operation ever
  decides key presence from a numerical zero-test.

This is (a) *mathematically identical* to C++ because the full operator is
block-diagonal across components, and (b) preserves the C++ complexity `Σ_g dim(g)³`,
mapping directly onto GPU batched solves later.

**Singular / near-singular policy.** The compatibility contract (equivalence within
tolerance) covers **non-singular inputs only**; all physical test fixtures are
well-conditioned, so this path is not exercised for backend-agreement. For a singular
component, behavior is *defined but not claimed bit-compatible with C++*: NumPy catches
`LinAlgError` and emits an **all-NaN component** (documented, backend-independent
semantics), with `xp.errstate` suppressing spurious warnings. This flows into the
existing NaN guard in `chiq_main.output2hdf5` exactly as a C++ inf/nan result would.
(Whether that guard should hard-error is a separate, out-of-scope fix.)

### 3.6 Files touched / added

- New: `python/package/chiq/solver/__init__.py` (factory), `base.py`, `layout.py`,
  `kernels.py`, `cpp.py`, `numpy.py`, `cupy.py`.
- Edit: `python/scripts/chiq_main.py` (use factory; drop direct `bse_solver` import
  and inline solver construction; drop `_reformat_data`/reshapes now in `layout`),
  `python/package/chiq/bse_toml.py` (parse+validate `backend`).
- Edit (packaging): `src/bse_solver_pybind.cpp` — rename the module declaration
  `PYBIND11_MODULE(bse_solver, m)` → `PYBIND11_MODULE(_bse_solver, m)` so the compiled
  binary exports `PyInit__bse_solver` and is importable as `chiq._bse_solver`. (Only
  the module name changes; the `BSESolver` class and its bindings are identical.)
  `src/CMakeLists.txt` / `python/CMakeLists.txt` install the extension into the `chiq`
  package directory. Both the pip (scikit-build-core) and direct-CMake paths install to
  the same package-relative location, so `import chiq._bse_solver` works identically
  (§5.3, verified by §4.5 tests). The top-level name `bse_solver` survives as a
  pure-Python shim (§5.2).

## 4. Testing / numerical verification

Layers:

1. **Backend-agreement (primary).** Parametrize the existing integration fixtures
   (`tests/python/non-mpi/bsetool_{BSE,RPA,SCL,RRPA}`) over `backend ∈ {cpp, numpy}`;
   assert `numpy` == `cpp` on identical input. Two fixed, non-discretionary metrics,
   aggregated as the max over all present blocks (a NaN/Inf in one backend but not the
   other, at the same position, is an immediate failure):
   - **mixed-tolerance absolute error** (standard `numpy.allclose` convention):
     `max|cpp - numpy| <= atol + rtol*max|cpp|`, with `atol=1e-10`, `rtol=1e-8`. This
     permits ~`atol` absolute error where the reference is zero/near-zero; the earlier
     `Δ/(scale+atol)` form was wrong (it forced ~`1e-18` at zero) and is replaced.
   - **relative backward error** of each internal solve,
     `||(I-M)x - b||_inf / (||b||_inf + atol) <= 1e-8`, evaluated at a
     **kernel/linalg-adapter test seam** — not through the public facade, so no debug
     state is added to `set/calc/get`.
   Both thresholds are fixed constants in the test module. A tolerance change requires
   an explicit spec amendment (not an ad-hoc relaxation). The backward-error metric is
   platform-robust where a raw `allclose` on explicit inverses may not be.
2. **Reference data.** Existing `.dat`/HDF5 comparisons also run through `numpy`.
3. **Independent analytical unit tests** (not just C++ as oracle): each formula and
   each index transform (`Sumfreq_A/B`, `Convert_A2B/B2A`, per-component inverse with
   fill-in) on small hand-checkable matrices, including **non-commuting blocks** to
   catch multiplication-order errors and **multi-block / missing-block** (sparse
   component) cases to catch index permutations and fill-in mistakes. All fixtures use
   the uniform-dimension invariant (§3.4); non-uniform info is tested only to assert it
   is **rejected**. Fixtures are constructed well-conditioned (bounded condition
   number) so agreement thresholds are meaningful.
4. **Config tests:** `backend` missing (→ `cpp`), wrong case (normalized), invalid
   value (clear error), and every existing sample TOML still loads unchanged.
5. **Packaging tests:** build sdist + wheel; `pip install` each into a clean venv
   *outside* the repo; smoke-test `import chiq`, `import chiq._bse_solver`,
   `import bse`, extension load, every console-script `--version`/`--help`, and
   packaged data. Run for each supported Python.

The CuPy skeleton is covered only by a factory test asserting the clear "not enabled"
error.

## 5. Packaging (pip / scikit-build-core)

### 5.1 Module layout (wheel + editable)

- The pybind11 extension is installed **inside** the package as `chiq._bse_solver`
  (deterministic in wheels; no reliance on a top-level module or on `PYTHONPATH`).
  `chiq.solver.cpp` does `from chiq import _bse_solver`.
- **Top-level `bse_solver` shim**: a tiny pure-Python `bse_solver.py`
  (`from chiq._bse_solver import *`) preserves `import bse_solver` for one release.
- **`bse` compatibility**: replace the filesystem *symlink* (unreliable in wheels)
  with a real shim package `bse/`. **Supported import forms** (the only ones that
  exist in the current tree, enumerated so the shim is testable): `import bse`,
  `import bse.<sub>` and `from bse import <sub>` for each existing submodule
  (`h5bse`, `bse_toml`, `matrix_dict`, `point_group`, `sumk_dft_chi`, `tools`,
  `index_pair`, `mpi`, …). Implemented as an explicit `bse/__init__.py` that, for each
  known submodule name, registers `sys.modules['bse.<sub>'] = importlib.import_module
  ('chiq.<sub>')` and re-exports them — deterministic, import-order-independent, no
  package-identity guessing. Emits one `DeprecationWarning` on first import; removed in
  the following release (removal version named in the deprecation message). Covered by
  module-identity and import-order tests (§4.5).

### 5.2 Entry points

- `console_scripts` for the 10 scripts that already expose `main()`:
  `chiq_main`, `chiq_post`, `chiq_fft`, `gen_qpath`, `gen_allq`, `calc_Iq`,
  `calc_Iq_scl`, `plot_chiq_path`, `plot_Ir`, `eigenvec_viewer`.
- `dcore_chiq.py` currently has only an `__main__` block; extract a `main()` and add
  its entry point.
- For one release, the legacy `<name>.py` command spellings are **also** shipped as
  aliases (documented as deprecated). `chiq_post` imports of `chiq_main` become a
  package-internal import rather than same-directory script import.

### 5.3 Build backend

- `pyproject.toml` with **scikit-build-core**: `pip install .` (and `-e .`) drives the
  existing top-level `CMakeLists.txt`, builds `bse_solver`, and lands it as
  `chiq._bse_solver`. The direct-CMake build path is kept working and **both paths are
  covered by the §4.5 install tests** so extension destination/RPATH/command set can't
  silently diverge (Codex risk).
- **Artifact policy (single source of truth).** The runtime **wheel** contains only
  genuine runtime resources: the `chiq` package, the `chiq._bse_solver` extension, and
  `point_group_data` (small, imported at runtime). Test `.h5` fixtures are **excluded
  from the wheel** and shipped only in the **sdist** (and referenced by the `test`
  extra / the test tree). Two distinct assertions: the **wheel** smoke test asserts
  `import chiq`, extension load, and `point_group_data` load (and does *not* look for
  test `.h5`); the **sdist**/test job asserts the `.h5` fixtures are present and load.

### 5.4 Dependency / Python matrix (pinned)

- `requires-python = ">=3.10"`. CI runs the **full build+ctest+pytest** on `3.10` and
  `3.13` (the range endpoints), plus a **lightweight wheel/import/backend smoke job**
  on `3.11` and `3.12`, so the whole advertised range is exercised.
- Build: `scikit-build-core>=0.10`, `pybind11>=2.12` (2.12+ supports 3.13),
  `cmake>=3.15`. Runtime deps are pinned by a **per-script import audit**: core =
  `numpy>=1.23`, `more-itertools`, `h5py` (HDF5 I/O in `h5bse`), `tomli` (only for
  <3.11; 3.11+ uses stdlib `tomllib`). Optional extras: `mpi[mpi4py]`,
  `dcore` (χ₀ preprocessing scripts only), `plot[matplotlib]` (plotting scripts only),
  `test[pytest]`, `gpu[cupy]` (later). Each advertised console script is smoke-tested
  in the extra that provides its deps.
- If a verified 3.13 wheel/source build of any dep is not yet available at
  implementation time, 3.13 drops to "best-effort" in CI and the floor/`3.10` target
  is unaffected.

### 5.5 MPI note

All ranks import `chiq._bse_solver` and construct the solver identically; `backend`
comes from the shared TOML so selection is uniform across ranks. No per-rank divergence.

## 6. Rollout / order of work

1. `chiq.solver` package: `layout` + `SolverBase` + `CppSolver`; route `chiq_main.py`
   through the factory (default `cpp`). Pure refactor — all tests stay green.
2. `kernels` + `NumpySolver`; backend-agreement + analytical tests; drive cpp↔numpy
   to tolerance.
3. `CupySolver` skeleton + factory wiring + "not enabled" test.
4. pip packaging (pyproject + scikit-build-core + entry points + `bse`/`bse_solver`
   shims); wheel/sdist install tests; update CI and docs (`README.md`,
   `doc/install.rst`).

Each step is independently testable and leaves the tree green.

## 7. Later phases (context only, not this spec)

- Phase 1: sparse-ir in χ₀ (`sumk_dft_chi.py`, `dcore_chiq.py`) to remove the
  fermionic-frequency-cutoff systematic error.
- Phase 2: enable/tune the CuPy backend on GPU (1 rank/GPU, memory-aware).
- Phase 3: IR/DLR compression of the BSE internal frequency sum itself.

## 8. Review resolution record (Codex rounds 1–3)

**Reviewer status.** Codex (technical/architecture lens) ran successfully both rounds.
Antigravity (user/workflow lens) is **non-functional in this environment**: the MCP
adapter invokes `agy -p --print-timeout 595s` without passing the design as a prompt,
so it returns an unrelated answer about the `--print-timeout` flag rather than a
review. Treated as an environment defect, not an approval; the user/workflow lens is
covered by the author + human review gate instead.

**Round-1 `must_fix_design` — all resolved:** groupwise-inverse invariant + per-
component solve (§3.5); state machine + per-calc contract (§3.3); canonical layout
owned by `chiq.solver.layout` (§3.4); ordered formulas, solve orientation, singular
policy (§3.5); wheel module layout, extension destination, `bse`/`bse_solver` shims,
entry points, install tests (§5).

**Round-2 `must_fix_design` — all resolved:** vertex-vs-key notation and per-conversion
key/shape contract (§3.4 table + reductions); uniform-dimension contradiction resolved
by committing to the uniform invariant and rejecting non-uniform, and dropping the
self-contradictory uneven-block test (§3.4, §4.3); structural fill-in rule for
inverse/solve/product matched to `block_matrix.hpp` (§3.5); executable singular-matrix
policy = all-NaN component, contract scoped to non-singular inputs (§3.5).

**Round-2 `should_fix` absorbed:** stale-output/no-caching policy (§3.3); mixed-layout
product contract via `ProdAB`/`Sumfreq` key equations (§3.4–3.5); fixed residual/error
norms and thresholds (§4.1); explicit enumerated `bse` shim import forms (§5.1);
`CMakeLists.txt` listed as edited for packaging (§3.6, §5.3); 3.11/3.12 smoke jobs and
per-script dependency audit (§5.4).

**Open questions — dispositions:**
- *matrix_info_* meaning* → resolved: per-vertex dimension arrays (§3.4).
- *inverse emits all intra-component pairs?* → resolved: yes, with numerical zeros
  retained (§3.5).
- *set-invalidates-outputs?* → accepted: no caching, `set→calc→get` ordering (§3.3).
- *key ordering / emitted zeros observable?* → resolved: keys are unordered; zeros
  retained within components; HDF5 write already prunes `|x|<1e-12`
  (`chiq_main.output2hdf5`), so on-disk output is unchanged (§3.5).
- *direct-CMake must keep top-level `bse_solver`?* → resolved: both build paths use the
  one canonical `chiq._bse_solver` layout; top-level name is a Python shim only (§5.1,
  §5.3).
- *`backend` TOML round-trip?* → resolved: new optional key in `[chiq_main]`, absent ⇒
  `cpp`; old parsers/generators that ignore unknown keys are unaffected (§3.2, §4.4).

**Risks — dispositions:**
- Ill-conditioned explicit inverses (Pq/I_q/Ueff/SCL) → *accepted*: inherent to the
  established C++ method; contract is equivalence on well-conditioned inputs, and
  backward-error metrics surface divergence (§4.1). Not introduced by this work.
- Sparsity-graph closure memory → *mitigated*: per-component solve caps at component
  size; components are small in ChiQ's block structure. *followup*: revisit if a future
  physical case yields one giant component (would also affect C++).
- Absent-key = exact zero assumption → *resolved*: §3.4 states absent ⇒ mathematically
  forbidden coupling; ChiQ builds dicts by structural rules, not by thresholding, so
  the assumption holds. Non-uniform/omitted-small-block inputs are rejected.
- C++ oracle hides shared conceptual errors → *mitigated*: independent analytical unit
  tests per formula/transform (§4.3).
- Test HDF5 in runtime wheel bloats size → *resolved*: §5.3 artifact policy — test
  `.h5` in sdist/test only, excluded from the runtime wheel; only `point_group_data`
  ships in the wheel; separate wheel vs sdist assertions.
- One-release deprecation breakage → *resolved*: removal version named in the
  `DeprecationWarning` for both `bse` and `bse_solver` shims (§5.1, §5.2).

**Round-3 `must_fix_design` — all resolved:**
- pybind module name → *resolved*: rename `PYBIND11_MODULE(bse_solver→_bse_solver)` so
  `PyInit__bse_solver` exists and `import chiq._bse_solver` works; top-level
  `bse_solver` becomes a pure-Python shim (§3.6, §5.1–5.2).
- tolerance formula error → *resolved*: switched to the standard mixed-tolerance rule
  `max|Δ| ≤ atol + rtol·max|cpp|` (§4.1).
- wheel-contents contradiction → *resolved*: single artifact policy with distinct
  wheel/sdist assertions (§5.3).
- product key-emission "numerical non-zero" → *resolved*: C++ `operator*` actually
  emits keys **symbolically** (`exists(i,k) && exists(k,j)`, `block_matrix.hpp:642`),
  not by value; §3.5 states the exact key-path rule, so the sparsity graph is
  backend-identical. Codex's premise (value-based) was mistaken; the corrected rule is
  strictly structural.

Round-3 `should_fix`/questions of note: `Convert_A2B/B2A` full index equations, exact
exception class for get-before-compute, negative/empty-dimension validation, and upper
dependency bounds are **deferred to TDD implementation** (they are concrete enough that
tests will pin them; carrying them as failing tests is more reliable than more prose).
The layout module's validation and the analytical unit tests (§4.3) are where these
become executable.
