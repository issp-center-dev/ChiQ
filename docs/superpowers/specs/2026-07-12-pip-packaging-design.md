# Design: pip packaging via scikit-build-core (coexisting with the CMake workflow)

**Date:** 2026-07-12
**Status:** Ready for review (round 2 — incorporates Codex + Antigravity design review)
**Repo:** ChiQ (Bethe-Salpeter / DMFT momentum-dependent susceptibility solver)
**Branch:** `design/pip-packaging` (stacked on `design/numpy-backend-packaging` / PR #6,
which adds the `chiq.solver` subpackage this design references)

## 1. Context and goal

ChiQ is currently installed **only** via CMake: `cmake && make && make install` copies the
`chiq` Python package to `<prefix>/lib/bse-python/`, installs the pybind11 extension
`bse_solver` next to it, drops the CLI scripts (`python/scripts/*.py`) into `bin`, and
generates `chiqvars.sh` which the user must `source` to set `PATH`/`PYTHONPATH`. There is
no `pyproject.toml`; `pip install` does not work.

**Goal:** make `pip install .` (and `pip install -e .`) build the extension and install the
package + CLI as a standard wheel, **without breaking the existing CMake + `chiqvars.sh`
workflow** — the two coexist. This is the packaging follow-up to the Phase 0 NumPy-backend
work; it distributes that work and lays groundwork for later optional deps (`gpu` → cupy,
IR → sparse-ir).

This design supersedes §5 of `2026-07-12-numpy-backend-and-packaging-design.md`: the Python
floor stays **3.8** (§6), and the CLI is restructured into an importable `chiq.cli`
subpackage (§4) so entry points are possible.

## 2. Non-goals / constraints

- **Coexistence, not replacement.** The direct `cmake && make install` + `chiqvars.sh`
  *workflow* (the commands a user runs) keeps working. `pip install .` is added alongside.
  Both are covered by install tests (§7). The two install *layouts* are selected by an
  explicit CMake option (`CHIQ_WHEEL_BUILD`, §3.1), not inferred from prefixes, so they
  cannot drift.
- **No behavior change** to the solver, I/O, or numerical CLI semantics. Packaging changes
  how the code is built, named on disk, imported, and installed — not what it computes.
- **Python floor stays 3.8** (ohtaka default; §6). Keep the third-party `toml` dep (stdlib
  `tomllib` is 3.11+) and 3.8-compatible syntax.
- **The C++ extension is pure Eigen — no MPI and no CUDA** (`grep` of `src/` finds no
  `mpi.h`/`MPI_`/`cuda`). MPI parallelism is Python-level (`mpi4py`, optional, runtime);
  GPU is the future Python CuPy backend. Therefore `pip` extras `mpi`/`gpu` are **runtime
  Python deps only** and never change C++ compile flags. No `-DUSE_MPI`/`-DUSE_GPU`.
- **Backward-compatible imports and commands** preserved, deprecated visibly, removed in a
  named future version (§10).

## 3. Build backend and native module layout

### 3.1 scikit-build-core config and the two install layouts

- Add `pyproject.toml` with **scikit-build-core** as the build backend, driving the existing
  top-level `CMakeLists.txt`.
- `[tool.scikit-build]` config:
  - `wheel.packages = ["python/package/chiq"]` — the `chiq` package is the wheel's importable
    tree (plus the shims, §3.3).
  - `cmake.args = ["-DCHIQ_WHEEL_BUILD=ON", "-DTesting=OFF"]` — a **new CMake option
    `CHIQ_WHEEL_BUILD`** (default `OFF`) selects wheel install destinations. When `ON`,
    `install(TARGETS _bse_solver)` targets the `chiq/` package dir (so the extension lands
    at `chiq/_bse_solver.<ext>`); the scripts/`chiqvars.sh`/`bse-python` install rules are
    skipped (the wheel gets its CLI from entry points, §4). When `OFF` (the legacy default),
    the current destinations are used **plus** the two new shims are installed into
    `lib/bse-python` (§3.3) so both layouts expose the same import surface.
  - `wheel.exclude`/package-data configured so the wheel contains only `chiq/` (incl.
    `chiq/point_group_data/` and `chiq/_bse_solver`), never `lib/bse-python/`, raw scripts,
    or `share/chiqvars.sh` (§7 artifact policy).
- Because the two layouts are chosen by one explicit option, `ordinary` CMake users are
  unaffected and the wheel layout is testable in isolation.

### 3.2 pybind11 source (submodule vs build dep)

The top-level `CMakeLists.txt` currently does `add_subdirectory(extern/pybind11)`
unconditionally. **An sdist does not contain the git submodule**, so a `pip install` from
sdist has no `extern/pybind11`. Resolution:

- Gate pybind11 acquisition on `CHIQ_WHEEL_BUILD`:
  - **Wheel/pip build** (`CHIQ_WHEEL_BUILD=ON`): `find_package(pybind11 CONFIG REQUIRED)`
    using the `pybind11` build-requirement installed into the PEP 517 build environment.
  - **Legacy direct build** (`OFF`): keep `add_subdirectory(extern/pybind11)` (the submodule).
- The pybind module symbol/target is renamed identically in both paths (§3.3), so the only
  difference is where pybind11's CMake config comes from.

### 3.3 Native module rename and backward-compat shims

- **Rename** `PYBIND11_MODULE(bse_solver, m)` → `PYBIND11_MODULE(_bse_solver, m)` in
  `src/bse_solver_pybind.cpp` (exports `PyInit__bse_solver`); set the CMake target
  `OUTPUT_NAME _bse_solver` so the `.so` basename matches. `BSESolver` and bindings unchanged.
  The extension installs as `chiq._bse_solver` in **both** layouts (wheel: `chiq/`; legacy:
  `lib/bse-python/chiq/`). `chiq/solver/cpp.py` does `from chiq import _bse_solver`.
- **Top-level `bse_solver` shim** (`python/package/bse_solver.py`,
  `from chiq._bse_solver import *`) preserves `import bse_solver; bse_solver.BSESolver`.
  Installed in **both** layouts (wheel: top-level module; legacy CMake: into `lib/bse-python`,
  a new install rule — this closes the gap that the legacy path previously installed only
  `chiq`). Emits a visible deprecation notice on import (§10).
- **`bse` package — lazy alias, not eager import** (`python/package/bse/__init__.py`):
  the previous "register every submodule in `sys.modules`" plan is unsound because it would
  eagerly import `chiq.sumk_dft_chi` (which imports `mpi4py` and `dcore`) and
  `chiq.g2scl_core` (matplotlib), breaking a core install. Instead the shim installs a
  **`MetaPathFinder`** that, on demand, resolves any `bse` or `bse.<dotted.path>` import to
  the corresponding `chiq.<dotted.path>` module and registers the *same module object* under
  the `bse.*` name (guaranteeing identity, incl. nested `bse.point_group_data.C1`). Nothing
  is imported until the user actually imports it, so optional deps are pulled only when the
  corresponding functionality is used. Supported forms (all tested, §7): `import bse`,
  `import bse.h5bse`, `from bse import h5bse`, `import bse.point_group_data.C1`,
  `importlib.import_module('bse.<name>')`, and mixed `chiq`/`bse` import order. Replaces the
  filesystem symlink (unreliable in wheels); the legacy CMake path installs this real shim
  package too (no more `bse -> chiq` symlink), so both layouts share one deprecation behavior.
  Emits one visible deprecation notice per process on first `bse` import.

## 4. CLI restructure into `chiq.cli` (required for entry points)

PEP 621 `console_scripts` must target **importable package modules**; the current
`python/scripts/*.py` are standalone files (not in any package) with same-directory sibling
imports (`chiq_post` → `import chiq_main`; `chiq_fft`, `gen_allq` → `from gen_qpath import
GenQPath`), which cannot be expressed as entry points. Resolution:

- **Move the CLI implementations into a new `chiq.cli` subpackage**: one module per command
  under `python/package/chiq/cli/` (`chiq_main.py`, `chiq_post.py`, `chiq_fft.py`,
  `gen_qpath.py`, `gen_allq.py`, `calc_Iq.py`, `calc_Iq_scl.py`, `plot_chiq_path.py`,
  `plot_Ir.py`, `eigenvec_viewer.py`, `dcore_chiq.py`), each exposing `def main()`.
  - `dcore_chiq.py` gains a `main()` (currently only an `if __name__` block).
  - Sibling imports become intra-package: `from chiq.cli import chiq_main`,
    `from chiq.cli.gen_qpath import GenQPath`. All three sibling-import sites are migrated.
- **Entry points** in `pyproject.toml`: `chiq_main = "chiq.cli.chiq_main:main"`, etc. — one
  per command. Optional-dependency commands (`dcore_chiq` needs `dcore`; `plot_chiq_path`,
  `plot_Ir` need `matplotlib`) still install as entry points but their heavy deps are
  imported **lazily inside `main()`** (or the module already imports a now-core dep — see §6),
  so `--help`/`--version` work on a core install and a missing extra yields a clear
  `"<command> requires the '<extra>' extra: pip install chiq[<extra>]"` error only when the
  functionality is actually invoked.
- **Legacy compatibility.** `python/scripts/<name>.py` are kept as **thin wrappers**
  (`from chiq.cli.<name> import main; main()`) so the direct-CMake path installs them to
  `bin` exactly as today. Additionally, deprecated `<name>.py` console-script aliases (e.g.
  `chiq_main.py = "chiq.cli._deprecated:chiq_main_py"`) are shipped for one release; each
  prints a deprecation line to `stderr` before delegating (so shell users see it, §10).

## 5. Optional-dependency handling (module-scope import audit)

Confirmed by grep of the actual source:

- `chiq.sumk_dft_chi`, `chiq.dcore_chiq_worker`, `chiq.cli.dcore_chiq` import **dcore**
  (and `mpi4py`) at module scope → gated behind the `dcore`/`mpi` extras; these are only
  reached by the χ₀ preprocessing command, imported lazily in that command's `main()`.
- `chiq.g2scl_core` imports **matplotlib** at module scope (`matplotlib.use('Agg')`). Because
  a non-optional core module imports it, **matplotlib is a core runtime dependency** (not a
  `plot`-only extra); the earlier "plot-only" classification was wrong. The `plot` extra is
  removed (matplotlib is core); plotting commands need nothing beyond core.
- `chiq._mpi` imports `mpi4py`, but `chiq.mpi` falls back to `chiq._no_mpi` when absent, so
  core commands (`chiq_main` etc.) run without `mpi4py`.

## 6. Dependencies and Python matrix (pinned)

- **Python:** `requires-python = ">=3.8"` (kept for ohtaka). CI runs the full
  build+test on **3.8 and 3.13**; a note states only these endpoints are continuously
  validated (a 3.11 wheel/import smoke job is added to catch mid-range breakage).
- **Core runtime deps:** `numpy>=1.23`, `scipy`, `more-itertools`, `h5py`, `toml`,
  **`matplotlib`** (imported by `chiq.g2scl_core`, §5).
- **Optional extras:** `mpi` → `mpi4py`; `dcore` → `dcore` (χ₀ preprocessing). `gpu` → `cupy`
  (future). (`plot` extra dropped — matplotlib is core; `test` → `pytest`.)
- **Build deps** (bounded so a future release cannot silently drop 3.8 from an isolated
  build): `scikit-build-core>=0.10,<1.0`, `pybind11>=2.12,<3.0`, `cmake>=3.15`. The sdist
  build is exercised under Python 3.8 in CI (§7) to catch build-dep drift early.

## 7. Testing

- **Wheel/sdist install tests (the crux of "pip works"):** build both an sdist and a wheel;
  `pip install` each into a **clean venv outside the repo** and smoke-test: `import chiq`,
  `import chiq._bse_solver` (extension loads), the `bse`/`bse_solver` shims incl. nested and
  mixed-order forms (§3.3), `point_group_data` load, and each **core** console command's
  `--version`/`--help` (`chiq_main`, `chiq_post`, `chiq_fft`, `gen_qpath`, `gen_allq`,
  `calc_Iq`, `calc_Iq_scl`, `plot_chiq_path`, `plot_Ir`, `eigenvec_viewer`). The
  extra-gated command `dcore_chiq --help` is tested only in a `[dcore]` env. Run under 3.8
  and 3.13.
- **sdist content assertion:** unpack the **sdist archive itself** (not the installed env,
  which builds a wheel) to assert test `.h5` fixtures are present in the sdist; assert the
  built **wheel** archive contains `chiq/` (with `_bse_solver` + `point_group_data`) and
  **no** test `.h5`, `lib/bse-python`, raw scripts, or `chiqvars.sh`.
- **Editable install:** `pip install -e .` then `import chiq._bse_solver` from a directory
  **outside** the repo (so in-tree `.so`/source can't mask a broken editable install);
  document that C++ source edits require an explicit `pip install -e . --no-build-isolation`
  rebuild (scikit-build-core does not auto-recompile on import).
- **Direct-CMake path unchanged:** a CI job still runs `cmake -DTesting=ON && make &&
  source chiqvars.sh && ctest -V && pytest`, proving the legacy workflow still works and
  installs the extension at the same `chiq._bse_solver` location.
- **Shim tests:** `bse`/`bse_solver` import-form/identity/import-order tests; a test that the
  deprecation notice fires once per process.

## 8. Files touched / added

- Add: `pyproject.toml`; `python/package/bse_solver.py` (shim); `python/package/bse/__init__.py`
  (lazy MetaPathFinder shim); `python/package/chiq/cli/` (one module per command +
  `_deprecated.py` wrappers).
- Edit: `src/bse_solver_pybind.cpp` (module rename); top-level `CMakeLists.txt` +
  `src/CMakeLists.txt` + `python/CMakeLists.txt` (add `CHIQ_WHEEL_BUILD`, pybind11 selection,
  install destinations/shims for both layouts); `python/package/chiq/solver/cpp.py`
  (`from chiq import _bse_solver`).
- Edit: `python/scripts/*.py` → thin wrappers delegating to `chiq.cli.*`.
- Edit: `.github/workflows/main.yml` (pip job + install/sdist/editable tests + 3.11 smoke);
  `README.md`, `doc/install.rst`.

## 9. HPC / offline build (ohtaka)

- Default `pip install .` uses PEP 517 **build isolation**, which downloads the build
  requirements (`scikit-build-core`, `pybind11`, `cmake`) from PyPI — this **fails on
  air-gapped cluster nodes**. `doc/install.rst` documents the offline path:
  `pip install scikit-build-core pybind11 cmake` (or `module load` them) into the target
  env, then `pip install . --no-build-isolation`.
- Compiler / Eigen selection: cluster users pass toolchain and Eigen hints via
  `CMAKE_ARGS="-DEIGEN3_DIR=... -DCMAKE_CXX_COMPILER=..."` env var (scikit-build-core forwards
  it). Documented alongside. (This is also why the CMake workflow is retained — it is the
  recommended path when heavy custom compile configuration is needed.)
- Install-precedence note: `doc/install.rst` warns against having both a pip-installed
  `chiq` and a `chiqvars.sh`-sourced `PYTHONPATH` entry active at once (pick one env), since
  a stale top-level `bse_solver` on `PYTHONPATH` could otherwise mask the shim.

## 10. Deprecation timeline (visible, named version)

- The `bse` shim, the top-level `bse_solver` shim, and the `<name>.py` console-script aliases
  are deprecated in this release and **removed in 2.0** (named explicitly, not "next minor").
  Research users run long-lived scripts, so the window is a full major version.
- **Visibility:** default Python filtering hides `DeprecationWarning` from non-`__main__`
  code, so the shims emit a **`FutureWarning`** (shown by default) *and* the `.py` command
  aliases print a one-line deprecation notice to `stderr` before delegating — a shell user
  running `chiq_main.py` will see it. Each message names the replacement (`chiq_main`,
  `import chiq...`) and the removal version (2.0).

## 11. Review resolution (Codex + Antigravity, round 1)

- **CLI not importable / sibling imports (Codex):** resolved by moving CLI into `chiq.cli`
  with intra-package imports and `chiq.cli.<name>:main` entry points; legacy `scripts/*.py`
  become thin wrappers (§4).
- **Eager `bse` shim pulls optional deps (Codex):** resolved by a lazy `MetaPathFinder` that
  imports nothing until used and preserves nested-module identity (§3.3).
- **scikit-build-core/CMake install contract + pybind11 selection + sdist-has-no-submodule
  (Codex):** resolved by the `CHIQ_WHEEL_BUILD` option, explicit wheel destinations/packages,
  and `find_package(pybind11)` for pip builds (§3.1–3.2).
- **Shims must exist in both layouts / bse symlink inconsistency (Codex):** legacy CMake now
  installs the same real shims (no symlink) (§3.3).
- **matplotlib is core, not plot-only (Codex):** reclassified; `plot` extra dropped (§5–6).
- **Optional-dep CLI `--help` can't pass on base install (Codex):** lazy imports in the
  extra-gated commands; core commands smoke-tested on base, `dcore_chiq` only under `[dcore]`
  (§4, §7).
- **3.8 not enforced by unbounded build reqs (Codex):** bounded build deps + sdist-build test
  on 3.8 (§6).
- **C++ MPI/GPU compile flags (Antigravity):** N/A — extension is pure Eigen; extras are
  runtime-only; clarified (§2).
- **PEP 517 build isolation fails offline on HPC (Antigravity):** documented
  `--no-build-isolation` + pre-install path for ohtaka (§9).
- **Silent DeprecationWarning + aggressive timeline (both):** FutureWarning + stderr notice,
  removal deferred to a named 2.0 (§10).
- Editable-install semantics, sdist-archive assertion, package-data for `point_group_data`,
  intermediate-Python note (Codex should-fix): folded into §6–§7.
