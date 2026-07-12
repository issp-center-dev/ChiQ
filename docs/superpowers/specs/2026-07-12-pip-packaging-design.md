# Design: pip packaging via scikit-build-core (coexisting with the CMake workflow)

**Date:** 2026-07-12
**Status:** Ready for review
**Repo:** ChiQ (Bethe-Salpeter / DMFT momentum-dependent susceptibility solver)
**Branch:** `design/pip-packaging` (stacked on `design/numpy-backend-packaging` / PR #6,
which adds the `chiq.solver` subpackage this design references)

## 1. Context and goal

ChiQ is currently installed **only** via CMake: `cmake && make && make install` copies
the `chiq` Python package to `<prefix>/lib/bse-python/`, installs the pybind11 extension
`bse_solver` next to it, drops the CLI scripts into `bin`, and generates `chiqvars.sh`
which the user must `source` to set `PATH`/`PYTHONPATH`. There is no `pyproject.toml`;
`pip install` does not work.

**Goal:** make `pip install .` (and `pip install -e .`) build the extension and install
the package + CLI as a standard wheel, **without breaking the existing CMake +
`chiqvars.sh` workflow** — the two coexist. This is the packaging follow-up to the Phase 0
NumPy-backend work; it distributes that work and lays groundwork for later optional deps
(`gpu` → cupy, IR → sparse-ir).

This design supersedes §5 of `2026-07-12-numpy-backend-and-packaging-design.md` on one
point: the Python floor stays at **3.8** (see §5), not 3.10.

## 2. Non-goals / constraints

- **Coexistence, not replacement.** The direct `cmake && make install` + `chiqvars.sh`
  *workflow* (the commands a user runs) keeps working. `pip install .` is added alongside.
  Both are covered by install tests (§6). Note the CMake install *layout* is adjusted so the
  extension lands at the same `chiq._bse_solver` location as the pip build (§9) — the
  commands and `chiqvars.sh` usage are unchanged, only where the `.so` is placed moves.
- **No behavior change** to the solver, I/O, or CLI semantics. Packaging only changes how
  the code is built, named on disk, and installed.
- **Python floor stays 3.8.** ohtaka (a production cluster) ships Python 3.8 by default;
  dropping it would break real users. Consequences: keep the third-party `toml` dependency
  (stdlib `tomllib` is 3.11+), and keep the Python code 3.8-compatible (no 3.9+ syntax).
- **Backward-compatible imports and commands** preserved for one release, then removed with
  a named-version `DeprecationWarning`.

## 3. Build backend and module layout

### 3.1 Build backend

- Add `pyproject.toml` with **scikit-build-core** as the build backend. `pip install .`
  and `pip install -e .` drive the existing top-level `CMakeLists.txt`, build the pybind11
  extension, and install `chiq` + the extension + the CLI entry points.
- `[build-system] requires = ["scikit-build-core>=0.10", "pybind11>=2.12", "cmake>=3.15"]`,
  `build-backend = "scikit_build_core.build"`.
- scikit-build-core is configured (via `[tool.scikit-build]`) to point at the existing
  `CMakeLists.txt`, install only the wheel-relevant targets, and pass through the CMake
  options needed for the extension. The **`Testing` CMake option stays OFF** for the pip
  build (the C++ gtest target is not needed in a wheel).

### 3.2 Native module layout (the one C++-touching change)

- The compiled extension is installed **inside the package as `chiq._bse_solver`**
  (deterministic in wheels; no reliance on a top-level module or `PYTHONPATH`). This
  requires:
  - `src/bse_solver_pybind.cpp`: rename `PYBIND11_MODULE(bse_solver, m)` →
    `PYBIND11_MODULE(_bse_solver, m)` (exports `PyInit__bse_solver`). The `BSESolver`
    class and all bindings are unchanged.
  - `src/CMakeLists.txt`: rename the pybind target (or set `OUTPUT_NAME _bse_solver`) so
    the `.so` basename is `_bse_solver.<ext>` and matches `PyInit__bse_solver`; install it
    into the `chiq` package directory.
  - `chiq/solver/cpp.py`: change `import bse_solver` → `from chiq import _bse_solver` and
    use `_bse_solver.BSESolver`.
- The leading underscore marks the raw binding as an internal implementation detail behind
  the `chiq.solver` facade (ecosystem convention: numpy `_core`, scipy `_fblas`, …). Users
  select a backend through `chiq.solver.get_solver` / the TOML `backend` option, never by
  importing the binding directly.

### 3.3 Backward-compatibility shims (one release)

- **Top-level `bse_solver`**: a tiny pure-Python `bse_solver.py`
  (`from chiq._bse_solver import *`) so the legacy `import bse_solver; bse_solver.BSESolver`
  keeps working. Emits a `DeprecationWarning` on import naming the removal version.
- **`bse` package**: today `python/package/bse` is a filesystem *symlink* to `chiq`
  (unreliable in wheels). Replace it with a real shim package `bse/` whose `__init__.py`
  registers `sys.modules['bse.<sub>'] = importlib.import_module('chiq.<sub>')` for each
  existing submodule (`h5bse`, `bse_toml`, `matrix_dict`, `point_group`, `sumk_dft_chi`,
  `tools`, `index_pair`, `mpi`, …) and re-exports them, so `import bse` and
  `from bse import <sub>` keep working. Emits one `DeprecationWarning` (removal version
  named). The CMake install path keeps creating its symlink as today (unchanged).

## 4. Entry points (CLI)

- Replace "install `python/scripts/*.py` into `bin`" with **`console_scripts`** for the ten
  scripts that already expose `main()`: `chiq_main`, `chiq_post`, `chiq_fft`, `gen_qpath`,
  `gen_allq`, `calc_Iq`, `calc_Iq_scl`, `plot_chiq_path`, `plot_Ir`, `eigenvec_viewer`.
- `python/scripts/dcore_chiq.py` currently has only an `if __name__ == "__main__"` block;
  extract a `main()` and add its entry point.
- For one release, the legacy `<name>.py` command spellings (e.g. `chiq_main.py`) are also
  shipped as `console_scripts` aliases pointing at the same `main()`, documented as
  deprecated. `chiq_post`'s same-directory `import chiq_main` becomes a package-internal
  import (the scripts must be importable as modules for entry points — verify none rely on
  being run as a file path).
- The direct-CMake path keeps installing the raw `.py` files to `bin` as today (unchanged),
  so both workflows provide the commands.

## 5. Dependencies and Python matrix (pinned)

Derived from an **actual per-file import audit** (not assumed):

- **Python:** `requires-python = ">=3.8"` (kept for ohtaka). CI matrix stays **3.8 and
  3.13** (the current endpoints).
- **Core runtime deps:** `numpy>=1.23`, `scipy` (imported by `chiq.g2scl_core` and the
  `chiq_fft`/`calc_Iq`/`calc_Iq_scl` scripts), `more-itertools`, `h5py` (HDF5 I/O in
  `h5bse`), and **`toml`** (the third-party package `bse_toml.py` imports for read *and*
  write — **not** `tomllib`, which is 3.11+ only; required to keep 3.8).
- **Optional extras** (PEP 621 strings pinned in `pyproject.toml`): `mpi` → `mpi4py`;
  `dcore` → `dcore` (χ₀ preprocessing scripts only); `plot` → `matplotlib` (plotting
  scripts only); `test` → `pytest`; `gpu` → `cupy` (reserved for the future GPU backend).
- **Build deps:** `scikit-build-core>=0.10`, `pybind11>=2.12` (works on 3.8–3.13),
  `cmake>=3.15`.
- A CI import-audit step re-derives the runtime dep list so it cannot silently drift.

## 6. Testing

- **Wheel/sdist install tests (the crux of "pip works"):** build both an sdist and a wheel;
  `pip install` each into a **clean virtual environment outside the repo** and smoke-test:
  `import chiq`, `import chiq._bse_solver` (extension loads), `import bse` and
  `from bse import h5bse` (shim), `import bse_solver` (top-level shim), each advertised
  console script's `--version`/`--help` (side-effect-free), and that `point_group_data`
  loads. Run per supported Python.
- **Direct-CMake path unchanged:** a job still runs `cmake -DTesting=ON && make &&
  source chiqvars.sh && ctest -V && pytest`, proving the legacy workflow still works and
  installs to the same package-relative extension location.
- **Module-identity / import-order tests** for the `bse` shim (each supported import form).
- **Deprecation-warning tests**: importing `bse` / `bse_solver` emits exactly one
  `DeprecationWarning` naming the removal version.
- **Artifact policy** (single source of truth): the runtime **wheel** contains only the
  `chiq` package, the `chiq._bse_solver` extension, and `point_group_data`. Test `.h5`
  fixtures are **excluded from the wheel**, shipped only in the **sdist**. Two distinct
  assertions: the wheel smoke test asserts no test `.h5` is present; the sdist/test job
  asserts they are.

## 7. CI and docs

- **CI (`.github/workflows/main.yml`):** keep the existing cmake→make→`source
  chiqvars.sh`→ctest→pytest job. **Add** a job that does `pip install .[test]` (+ the
  system deps eigen/hdf5/openmpi as now) and runs `pytest`, plus the wheel/sdist
  clean-venv install smoke tests. Python 3.8 and 3.13.
- **Docs:** add a "pip install" path to `doc/install.rst` and `README.md` (e.g.
  `pip install .` for a normal install, `pip install -e .` for development, extras like
  `pip install .[mpi,dcore]`), **keeping** the existing CMake instructions alongside. Note
  that after a pip install the console commands are on `PATH` without sourcing
  `chiqvars.sh`.

## 8. Files touched / added

- Add: `pyproject.toml`.
- Add: `python/package/bse_solver.py` (top-level shim), `python/package/bse/__init__.py`
  (shim package; replaces the symlink for the pip layout — the CMake symlink path stays).
- Edit: `src/bse_solver_pybind.cpp` (module rename), `src/CMakeLists.txt` (target
  name/`OUTPUT_NAME` + install destination for the pip build), `python/CMakeLists.txt` if
  needed for the coexisting install.
- Edit: `python/package/chiq/solver/cpp.py` (`from chiq import _bse_solver`; keep working
  under both install layouts — see §9).
- Edit: `python/scripts/dcore_chiq.py` (extract `main()`).
- Edit: `.github/workflows/main.yml` (add pip job + install tests), `README.md`,
  `doc/install.rst`.
- Add: install-test scripts/fixtures under `tests/` (clean-venv wheel/sdist smoke).

## 9. Key risk: `cpp.py` import under both layouts

`chiq/solver/cpp.py` must import the extension whether it was installed by pip (as
`chiq._bse_solver`) or by the legacy CMake path (which historically put a top-level
`bse_solver` on `PYTHONPATH`). Resolution: the CMake install is updated to also place the
extension as `chiq._bse_solver` (same package-relative location as pip), and `cpp.py` does
`from chiq import _bse_solver`. The top-level `bse_solver` name is provided by the
pure-Python shim in both layouts. This keeps a single import path in `cpp.py` and one
canonical on-disk location, avoiding a try/except-import fork.

## 10. Deprecation timeline

The `bse` package shim, the top-level `bse_solver` shim, and the `<name>.py` console-script
aliases are all deprecated in this release and **removed in the next minor release**; the
exact removal version is named in each `DeprecationWarning` message. (Version string lives
in `python/package/chiq/__init__.py`; currently `1.0-beta`.)
