# Codex handoff — finish ChiQ pip packaging (Plan B)

You are taking over work on **ChiQ** (a C++/pybind11 + Python scientific code: Bethe–Salpeter
equation solver within DMFT). Prior work was done by Claude Code. Your job is to **execute
"Plan B" (build & packaging)** to make `pip install .` work, then open a PR. Everything you
need is in the repo; this file is your orientation.

Work in: `/Users/k-yoshimi/Dropbox/cursor/ChiQ`. Current branch: **`design/pip-packaging`**
(do not switch branches).

---

## 1. What is already done (do NOT redo)

Two stacked branches off `main`:

- **`design/numpy-backend-packaging`** → GitHub **PR #6** (open): added a
  `backend = "cpp"|"numpy"|"cupy"` switch (the `chiq.solver` subpackage). Reviewed & green.
- **`design/pip-packaging`** (current, stacked on the above; **not pushed**): the pip-packaging
  effort, split into two plans:
  - **Plan A — DONE** (Python restructure, fully reviewed, 101/101 non-mpi tests green):
    - `chiq.g2scl_core` imports matplotlib lazily.
    - `python/package/bse` symlink → real **static forwarding-shim package** (`bse/<mod>.py`
      = `from chiq.<mod> import *`, allowlisted, lazy, emits `FutureWarning`).
    - CLI moved into **`chiq.cli`** package (11 command modules, each `def main()`);
      `dcore_chiq` has a lazy `_import_dcore()` + deferred class factory; plot/dcore commands
      import optional deps lazily inside `main()`.
    - `python/scripts/<name>.py` are now thin wrappers → `chiq.cli.<name>:main`;
      `python/package/chiq/cli/_deprecated.py` has `<name>_py` alias callables.
  - **Plan B — YOUR JOB, not started.**

## 2. Your task: execute Plan B

**Read the plan and do it task by task:**
`docs/superpowers/plans/2026-07-12-pip-packaging-B-build.md` (5 tasks, TDD, each with exact
file contents/commands and a commit).
**The approved design it implements:** `docs/superpowers/specs/2026-07-12-pip-packaging-design.md`.

Plan B summary (details in the plan):
1. **B1** Rename `PYBIND11_MODULE(bse_solver)` → `_bse_solver`; install it as `chiq._bse_solver`;
   update `chiq/solver/cpp.py` to `from chiq import _bse_solver`; add top-level `bse_solver/`
   shim package.
2. **B2** Add `CHIQ_WHEEL_BUILD` CMake option; pip build uses `find_package(pybind11)` (sdists
   have no git submodule), legacy build keeps `add_subdirectory(extern/pybind11)`; dual install
   layouts; **remove the old `bse` symlink** and install the real `bse`/`bse_solver` shims in the
   legacy layout.
3. **B3** `pyproject.toml` (scikit-build-core, deps, `[project.scripts]` entry points incl. the
   deprecated `"<name>.py"` aliases, extras).
4. **B4** Clean-venv **wheel / sdist / editable** install verification script.
5. **B5** CI job (`pip install .[test,plot]` on 3.8 + 3.13) + docs (README, `doc/install.rst`,
   `doc/tutorial/*.rst`).

Follow the plan's TDD steps and commit after each task. Keep commits scoped as the plan shows.

## 3. CRITICAL: local build environment (Apple clang 21 + Homebrew Eigen)

This machine's toolchain needs specific flags or the C++ build fails. **Export these before any
`pip install` / `cmake` / `make` in this repo** (scikit-build-core forwards `CMAKE_ARGS`):

```bash
export SDK=/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk
export CMAKE_ARGS="-DCMAKE_POLICY_VERSION_MINIMUM=3.5 \
  -DEIGEN3_INCLUDE_DIR=$(brew --prefix eigen)/include/eigen3 \
  -DPYTHON_INCLUDE_DIR=/Users/k-yoshimi/miniconda3/include/python3.10 \
  -DPYTHON_LIBRARY=/Users/k-yoshimi/miniconda3/lib/libpython3.10.a \
  -DCMAKE_CXX_FLAGS=-isysroot\ $SDK\ -isystem\ $SDK/usr/include/c++/v1"
```

Why each flag:
- `CMAKE_POLICY_VERSION_MINIMUM=3.5` — CMake here is 4.x, which rejects the project's old
  `cmake_minimum_required(3.5)` without this.
- `EIGEN3_INCLUDE_DIR` — Eigen was `brew install eigen`'d; it is not on the default search path.
- `PYTHON_INCLUDE_DIR`/`PYTHON_LIBRARY` — the bundled/older pybind11 FindPythonLibs otherwise
  picks a non-existent `python.app` include path.
- `-isysroot ... -isystem .../c++/v1` — **Apple clang 21 (clang-2100) has a broken default C++
  header search**; without these, `#include <iostream>` fails ("file not found").

Also note the project already requires **C++17** (`CMAKE_CXX_STANDARD 17` in `CMakeLists.txt`,
set in earlier work) because modern Eigen 3.4 needs ≥ C++14. Do not revert that.

**These flags are LOCAL ONLY — never commit them into `pyproject.toml`.** On Linux CI they are
unnecessary (Eigen via `libeigen3-dev`, normal clang/gcc, `find_package(pybind11)` from the
build-dep). `pyproject.toml`'s `cmake.args` should contain only portable flags
(`-DCHIQ_WHEEL_BUILD=ON -DTesting=OFF -DCMAKE_POLICY_VERSION_MINIMUM=3.5`).

A legacy build already exists at `build/`. The built extension after B1 lives at
`python/package/chiq/_bse_solver*.so` (gitignored). If the source-tree `PYTHONPATH` workflow
misbehaves after the rename, prefer `pip install -e .` for local test runs.

## 4. Test / verify commands

- **Python unit/integration tests** (before B1, extension is top-level `bse_solver` in `build/src`):
  ```bash
  export PYTHONPATH="$PWD/python/package:$PWD/build/src"
  export PATH="$PWD/python/scripts:$PATH"     # integration tests shell out to chiq_main.py etc.
  chmod +x python/scripts/*.py
  cd build && python3 -m pytest ../tests/python/non-mpi -q     # expect all green (101 before Plan B)
  ```
  After B1 (extension is `chiq._bse_solver` inside the package), drop `build/src` from
  `PYTHONPATH` (use `PYTHONPATH="$PWD/python/package"` + the copied `.so`, or an editable install).
- **Packaging unit tests:** `python3 -m pytest tests/python/non-mpi/packaging -q`.
- **pyflakes** is authoritative for "undefined name" bugs (a real one — a missing
  `from itertools import product` — was caught in Plan A's `dcore_chiq.py`; keep it clean):
  `python3 -m pyflakes python/package/chiq/cli/dcore_chiq.py` → no "undefined name".
- **Install verification** (Plan B B2/B4 create these): `bash
  tests/python/non-mpi/packaging/verify_legacy_install.sh` and `.../verify_pip_install.sh`.

## 5. Non-negotiable constraints

- **Python floor 3.8** (the ohtaka cluster's default). `requires-python=">=3.8"`, CI 3.8+3.13,
  no 3.9+ syntax, **keep the third-party `toml`** dep (it reads *and writes* TOML; stdlib
  `tomllib` is 3.11+ and read-only).
- **Coexistence:** the existing `cmake && make install` + `chiqvars.sh` workflow must keep
  working. Layout is selected by the `CHIQ_WHEEL_BUILD` CMake option (default OFF = legacy).
- **The C++ extension is pure Eigen — no MPI/CUDA.** `mpi`/`gpu` extras are runtime Python deps
  only; do NOT add `-DUSE_MPI`/`-DUSE_GPU`.
- **Extension = `chiq._bse_solver`** in both layouts; top-level `bse_solver` is a pure-Python
  shim. `matplotlib` is a `plot` extra (lazy). Deprecation via `FutureWarning`/stderr, removal
  in **2.0**.
- Don't change what the CLI commands compute — packaging only.

## 6. Review & finish

- The prior workflow ran per-task code reviews and a final review. You can self-review, or use
  the `ai-review-cycle` MCP tools if available in your Codex setup
  (`codex_review_diff` / `antigravity_review_diff` — **both reviewers work**; `agy`/Antigravity
  was broken earlier but is fixed as of 2026-07-12, verify with `doctor_reviewers`). Whatever you
  use, verify: `pip install .` builds and imports on a clean venv; the legacy `make install` +
  `chiqvars.sh` path still works; the non-mpi suite stays green.
- **Open the PR when done:** the branch is stacked on PR #6, so target `main` and note the stack
  (or wait until #6 merges, then rebase onto main). **`gh` needs `env -u GITHUB_TOKEN`** here
  (a stale `GITHUB_TOKEN` env var shadows the valid keyring auth):
  ```bash
  git push -u origin design/pip-packaging
  env -u GITHUB_TOKEN gh pr create --base main --head design/pip-packaging --title "..." --body "..."
  ```

## 7. Known gotchas / deferred items already captured in Plan B

- `python/CMakeLists.txt` still creates a `bse -> chiq` **symlink** at install; B2 removes it and
  installs the real shim (Antigravity flagged the symlink would shadow the static shim).
- Docs still use deprecated `.py` command names (`chiq_main.py`, etc.) and README has a stale
  `../bse/doc` sphinx path and a `chi_main.py` typo in `single_orb.rst`; B5 fixes all of these
  (full inventory is in spec §7-§8).
- If `pip install .` genuinely cannot build here despite the `CMAKE_ARGS`, that's a
  toolchain blocker, not a plan defect: capture the exact error; the legacy `make install`
  verification (B2) independently proves the packaging layout.

Good luck. Start by reading `docs/superpowers/plans/2026-07-12-pip-packaging-B-build.md` and doing
Task B1.
