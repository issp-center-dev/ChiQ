# AGENTS.md

Guidance for AI coding agents (Codex etc.) working in this repository. Keep it in sync with
`CLAUDE.md` (same project; that file has more architecture detail).

## What ChiQ is

C++/pybind11 + Python scientific code: solves the Bethe–Salpeter equation within DMFT and
computes momentum-dependent susceptibilities χ(q, iν_m). Runs as a post-processing step of
DCore. Two layers: a C++ Eigen solver exposed via pybind11 (`src/`), and a Python package
(`python/package/chiq/`) with CLI scripts.

## Build (this machine needs specific flags — Apple clang 21 + Homebrew Eigen)

Modern Eigen needs C++17 (already set in `CMakeLists.txt`). The local toolchain requires these
flags or the C++ build fails; **export before any `cmake`/`make`/`pip install`** (scikit-build-core
forwards `CMAKE_ARGS`). They are **local-only — never commit them into `pyproject.toml`**:

```bash
export SDK=/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk
export CMAKE_ARGS="-DCMAKE_POLICY_VERSION_MINIMUM=3.5 \
  -DEIGEN3_INCLUDE_DIR=$(brew --prefix eigen)/include/eigen3 \
  -DPYTHON_INCLUDE_DIR=/Users/k-yoshimi/miniconda3/include/python3.10 \
  -DPYTHON_LIBRARY=/Users/k-yoshimi/miniconda3/lib/libpython3.10.a \
  -DCMAKE_CXX_FLAGS=-isysroot\ $SDK\ -isystem\ $SDK/usr/include/c++/v1"
```
- `CMAKE_POLICY_VERSION_MINIMUM=3.5`: CMake 4.x rejects the project's old `cmake_minimum_required` without it.
- `EIGEN3_INCLUDE_DIR`: Eigen was `brew install`'d, not on the default path.
- `PYTHON_INCLUDE_DIR`/`PYTHON_LIBRARY`: pybind11 FindPythonLibs otherwise picks a non-existent `python.app` path.
- `-isysroot ... -isystem .../c++/v1`: Apple clang 21 has broken default C++ header search (`<iostream>` not found).

A legacy build tree exists at `build/`.

## Tests

```bash
export PYTHONPATH="$PWD/python/package:$PWD/build/src"   # after the _bse_solver rename, use just python/package (+ editable install)
export PATH="$PWD/python/scripts:$PATH"                  # integration tests shell out to chiq_main.py etc.
chmod +x python/scripts/*.py
cd build && python3 -m pytest ../tests/python/non-mpi -q # expect all green
```
- Packaging unit tests: `python3 -m pytest tests/python/non-mpi/packaging -q`.
- `python3 -m pyflakes <file>` is authoritative for "undefined name" bugs — keep changed files clean.
- Legacy C++ tests: `cmake -DTesting=ON && make && source chiqvars.sh && ctest -V`.

## Non-negotiable constraints

- **Python floor 3.8** (ohtaka cluster default): `requires-python=">=3.8"`, CI 3.8+3.13, no 3.9+
  syntax, **keep the `toml` dependency** (reads *and writes* TOML; stdlib `tomllib` is 3.11+, read-only).
- **Coexistence:** `pip install .` and the legacy `cmake && make install` + `chiqvars.sh` workflow
  must both work; the layout is selected by the `CHIQ_WHEEL_BUILD` CMake option (default OFF = legacy).
- **The C++ extension is pure Eigen — no MPI/CUDA compile flags.** `mpi`/`gpu` pip extras are
  runtime Python deps only.
- **Default behavior unchanged.** Solver backend switch defaults to `cpp`; packaging must not change
  what commands compute.
- The compiled extension is `chiq._bse_solver`; top-level `bse_solver` and the `bse` package are
  deprecated pure-Python shims (removed in ChiQ 2.0; warn via `FutureWarning`/stderr). `matplotlib`
  is a `plot` extra (imported lazily).

## Current work state (2026-07-12)

Branches off `main`:
- `design/numpy-backend-packaging` → **PR #6** (open): `backend="cpp"|"numpy"|"cupy"` solver switch.
- `design/pip-packaging` (stacked on #6, **not pushed**): pip packaging.
  - **Plan A DONE** (Python restructure: `chiq.cli` package, static `bse` shim, lazy matplotlib,
    script wrappers). Reviewed, non-mpi 101/101 green.
  - **Plan B TODO** (build & packaging): execute
    `docs/superpowers/plans/2026-07-12-pip-packaging-B-build.md`.

**If continuing the pip-packaging work, read `docs/superpowers/CODEX-HANDOFF-pip-packaging.md` first**
(full orientation), then execute Plan B task by task (TDD, commit per task). Design spec:
`docs/superpowers/specs/2026-07-12-pip-packaging-design.md`.

## Conventions

- Version string: `python/package/chiq/__init__.py` (`__version__`).
- `gh` needs `env -u GITHUB_TOKEN gh ...` here (a stale `GITHUB_TOKEN` env var shadows the valid
  keyring auth).
- Design docs, plans, and per-task ledgers live under `docs/superpowers/` and `.superpowers/sdd/`.
- Optional: the `ai-review-cycle` MCP provides Codex/Antigravity review tools (`codex_review_diff`,
  `antigravity_review_diff`, `doctor_reviewers`); both reviewers are working as of 2026-07-12.
