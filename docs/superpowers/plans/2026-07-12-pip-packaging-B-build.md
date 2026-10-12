# pip packaging — Plan B: build & packaging Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make `pip install .` build the pybind11 extension and install ChiQ as a standard wheel (with CLI entry points and backward-compat shims), coexisting with the existing `cmake && make install` + `chiqvars.sh` workflow. Depends on Plan A (`chiq.cli`, static `bse` shim).

**Architecture:** Add `pyproject.toml` (scikit-build-core) driving the existing `CMakeLists.txt`. Rename the pybind module to `chiq._bse_solver` and install it *inside* the `chiq` package. A new `CHIQ_WHEEL_BUILD` CMake option selects wheel vs legacy install destinations; pip builds use `find_package(pybind11)` (sdists have no git submodule), legacy builds keep the bundled submodule. Clean-venv wheel/sdist/editable install tests prove both paths.

**Tech Stack:** scikit-build-core, pybind11, CMake, Eigen3, pytest, Python ≥3.8.

## Global Constraints

- **Python floor 3.8** (ohtaka default). `requires-python = ">=3.8"`; CI 3.8 + 3.13; keep `toml` (read+write; `tomllib` is 3.11+ read-only). No 3.9+ syntax.
- **Coexistence:** the direct `cmake && make install` + `chiqvars.sh` workflow keeps working; `pip install .` is added. Layouts are chosen by the `CHIQ_WHEEL_BUILD` CMake option (default `OFF` = legacy), never inferred.
- **The C++ extension is pure Eigen — no MPI/CUDA compile flags.** `mpi`/`gpu` extras are runtime Python deps only.
- **Extension name:** `PYBIND11_MODULE(_bse_solver)`, installed as `chiq._bse_solver` in BOTH layouts. Top-level `bse_solver` survives as a pure-Python shim package. `chiq/solver/cpp.py` imports `from chiq import _bse_solver`.
- **matplotlib is a `plot` extra** (Plan A made it lazy). Core deps: `numpy>=1.23, scipy, more-itertools, h5py, toml`. Extras: `plot`→matplotlib, `mpi`→mpi4py, `dcore`→dcore, `gpu`→cupy, `test`→pytest.
- **Deprecation:** shims warn (`FutureWarning`/stderr), removed in **2.0** (named).
- Spec: `docs/superpowers/specs/2026-07-12-pip-packaging-design.md`.
- **LOCAL BUILD ENV (this dev machine — Apple clang 21 + Homebrew Eigen).** The pip build must be given the same flags the direct build needed. Export before any `pip install`/`make` in this repo:
  ```bash
  export SDK=/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk
  export CMAKE_ARGS="-DCMAKE_POLICY_VERSION_MINIMUM=3.5 \
    -DEIGEN3_INCLUDE_DIR=$(brew --prefix eigen)/include/eigen3 \
    -DPYTHON_INCLUDE_DIR=/Users/k-yoshimi/miniconda3/include/python3.10 \
    -DPYTHON_LIBRARY=/Users/k-yoshimi/miniconda3/lib/libpython3.10.a \
    -DCMAKE_CXX_FLAGS=-isysroot\ $SDK\ -isystem\ $SDK/usr/include/c++/v1"
  ```
  (scikit-build-core forwards `CMAKE_ARGS`. On a normal machine these are unnecessary — they are a local-toolchain workaround, NOT baked into `pyproject.toml`.) The already-built legacy tree is at `build/`.

## File Structure

- Modify: `src/bse_solver_pybind.cpp` (module rename); `src/CMakeLists.txt` (target rename + wheel/legacy install dest); top-level `CMakeLists.txt` (`CHIQ_WHEEL_BUILD` option, pybind11 selection, skip chiqvars in wheel); `python/CMakeLists.txt` (wheel skips; legacy installs `chiq` + `bse` + `bse_solver` shims, removes old symlink); `python/package/chiq/solver/cpp.py` (`from chiq import _bse_solver`).
- Create: `python/package/bse_solver/__init__.py` (top-level shim package); `pyproject.toml`.
- Create tests: `tests/python/non-mpi/packaging/test_install_wheel.py` (or a shell-driven CI job — see Task B4).
- Modify: `.github/workflows/main.yml`; docs (`README.md`, `doc/install.rst`, `doc/tutorial/*.rst`) per spec §7-§8.

---

### Task B1: Rename the pybind module to `_bse_solver` and add the top-level `bse_solver` shim

**Files:**
- Modify: `src/bse_solver_pybind.cpp` (the `PYBIND11_MODULE` line, ~line 223)
- Modify: `src/CMakeLists.txt` (target name + install)
- Modify: `python/package/chiq/solver/cpp.py`
- Create: `python/package/bse_solver/__init__.py`
- Test: `tests/python/non-mpi/packaging/test_bse_solver_shim.py`

**Interfaces:**
- Consumes: nothing new.
- Produces: extension importable as `chiq._bse_solver` (with `BSESolver`); top-level `import bse_solver` shim.

- [ ] **Step 1: Write the failing test**

```python
# tests/python/non-mpi/packaging/test_bse_solver_shim.py
import importlib
import sys
import warnings
import pytest

def test_chiq_has_underscore_bse_solver():
    m = importlib.import_module("chiq._bse_solver")
    assert hasattr(m, "BSESolver")

def test_top_level_bse_solver_shim_warns_and_forwards():
    sys.modules.pop("bse_solver", None)
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        import bse_solver
    assert hasattr(bse_solver, "BSESolver")
    assert any(issubclass(x.category, FutureWarning) for x in w)
```

- [ ] **Step 2: Run test to verify it fails**

Requires a rebuild with the extension placed as `chiq/_bse_solver`. First do Steps 3-4, then rebuild (Step 5) — this test cannot pass until the rebuilt `.so` is importable as `chiq._bse_solver`. Run after Step 5:
`PYTHONPATH="$(pwd)/python/package:$(pwd)/build/src" python3 -m pytest tests/python/non-mpi/packaging/test_bse_solver_shim.py -v` → initially FAIL (module named `bse_solver`, no `chiq._bse_solver`, no shim).

- [ ] **Step 3: Rename the pybind module**

In `src/bse_solver_pybind.cpp`, change (near line 223):
```cpp
PYBIND11_MODULE(bse_solver, m) {
```
to:
```cpp
PYBIND11_MODULE(_bse_solver, m) {
```
(Nothing else in the file changes — `BSESolver` and its bindings are identical.)

In `src/CMakeLists.txt`, change:
```cmake
pybind11_add_module(bse_solver bse_solver_pybind.cpp)

install(TARGETS bse_solver
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}/bse-python)
```
to:
```cmake
pybind11_add_module(_bse_solver bse_solver_pybind.cpp)

# Install the extension INSIDE the chiq package in both layouts (Task B2 sets
# CHIQ_WHEEL_BUILD; here we just target the chiq package dir).
if(CHIQ_WHEEL_BUILD)
  install(TARGETS _bse_solver LIBRARY DESTINATION chiq)
else()
  install(TARGETS _bse_solver
      LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}/bse-python/chiq)
endif()
```
(The `CHIQ_WHEEL_BUILD` option is defined in Task B2; until then it is undefined and the `else()` legacy branch is taken — that is correct for this task.)

- [ ] **Step 4: Update `cpp.py` and create the top-level shim**

`python/package/chiq/solver/cpp.py`: change
```python
import bse_solver  # top-level pybind module (packaging phase renames to chiq._bse_solver)
```
to
```python
from chiq import _bse_solver
```
and change the constructor's `bse_solver.BSESolver(...)` to `_bse_solver.BSESolver(...)`.

Create `python/package/bse_solver/__init__.py`:
```python
"""Backward-compatibility shim for the old top-level `bse_solver` extension name.

The compiled module is now `chiq._bse_solver`. This shim re-exports it so the
legacy `import bse_solver; bse_solver.BSESolver` keeps working. Deprecated;
removed in ChiQ 2.0 -- use `chiq.solver.get_solver` instead.
"""
import warnings

from chiq._bse_solver import *  # noqa: F401,F403

warnings.warn(
    "Importing the top-level 'bse_solver' module is deprecated and will be "
    "removed in ChiQ 2.0; use chiq.solver.get_solver(...) instead.",
    FutureWarning,
    stacklevel=2,
)
```

- [ ] **Step 5: Rebuild so the extension is importable as `chiq._bse_solver`, then run tests**

For the SOURCE-tree test workflow, the built `_bse_solver.<ext>` must sit inside `python/package/chiq/`. Rebuild and copy it there (the definitive test is the install tests in Task B4; this copy makes the source-tree PYTHONPATH workflow work for the unit tests):
```bash
export SDK=/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk
rm -rf build && mkdir build && cd build
cmake .. -DTesting=OFF -DCMAKE_BUILD_TYPE=Release -DCMAKE_POLICY_VERSION_MINIMUM=3.5 \
  -DEIGEN3_INCLUDE_DIR="$(brew --prefix eigen)/include/eigen3" \
  -DPYTHON_EXECUTABLE="$(which python3)" \
  -DPYTHON_INCLUDE_DIR=/Users/k-yoshimi/miniconda3/include/python3.10 \
  -DPYTHON_LIBRARY=/Users/k-yoshimi/miniconda3/lib/libpython3.10.a \
  -DCMAKE_CXX_FLAGS="-isysroot $SDK -isystem $SDK/usr/include/c++/v1"
make _bse_solver
cd ..
cp build/src/_bse_solver.*.so python/package/chiq/    # so `from chiq import _bse_solver` works from the source tree
```
Then run:
```bash
PYTHONPATH="$(pwd)/python/package" python3 -m pytest tests/python/non-mpi/packaging/test_bse_solver_shim.py -v
```
Expected: PASS (2 tests).
Also confirm the cpp backend still works via the factory:
```bash
PYTHONPATH="$(pwd)/python/package" python3 -c "from chiq.solver import get_solver; s=get_solver('cpp',1.0,[1],[1],[1]); print('cpp backend OK')"
```
And run the full non-mpi suite with the new layout (note: the built `.so` now lives in `python/package/chiq/`, so PYTHONPATH no longer needs `build/src`):
```bash
chmod +x python/scripts/*.py
export PYTHONPATH="$(pwd)/python/package"; export PATH="$(pwd)/python/scripts:$PATH"
cd build && python3 -m pytest ../tests/python/non-mpi -q   # all green
```
(`python/package/chiq/*.so` is gitignored — confirm `.gitignore` covers `*.so`; if not, add `python/package/chiq/_bse_solver*.so` to `.gitignore` in this commit.)

- [ ] **Step 6: Commit**

```bash
git add src/bse_solver_pybind.cpp src/CMakeLists.txt python/package/chiq/solver/cpp.py python/package/bse_solver/__init__.py tests/python/non-mpi/packaging/test_bse_solver_shim.py .gitignore
git commit -m "feat(build): rename pybind module to chiq._bse_solver + top-level bse_solver shim"
```

---

### Task B2: `CHIQ_WHEEL_BUILD` option, pybind11 selection, dual install layouts

**Files:**
- Modify: top-level `CMakeLists.txt`
- Modify: `python/CMakeLists.txt`
- Test: a legacy `make install` into a temp prefix (shell verification below)

**Interfaces:**
- Consumes: the `_bse_solver` target (B1).
- Produces: `CHIQ_WHEEL_BUILD` CMake option; wheel mode installs only the extension to `chiq/`; legacy mode installs `chiq` + `bse` + `bse_solver` shims into `lib/bse-python` and removes the old `bse` symlink.

- [ ] **Step 1: Write the failing test (a shell verification script)**

Create `tests/python/non-mpi/packaging/verify_legacy_install.sh`:
```bash
#!/usr/bin/env bash
# Verifies the legacy (CHIQ_WHEEL_BUILD=OFF) install layout into a temp prefix.
set -euo pipefail
PREFIX=$(mktemp -d)
SDK=/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk
rm -rf build_legacy && mkdir build_legacy && cd build_legacy
cmake .. -DTesting=OFF -DCMAKE_BUILD_TYPE=Release -DCMAKE_POLICY_VERSION_MINIMUM=3.5 \
  -DCMAKE_INSTALL_PREFIX="$PREFIX" \
  -DEIGEN3_INCLUDE_DIR="$(brew --prefix eigen)/include/eigen3" \
  -DPYTHON_EXECUTABLE="$(which python3)" \
  -DPYTHON_INCLUDE_DIR=/Users/k-yoshimi/miniconda3/include/python3.10 \
  -DPYTHON_LIBRARY=/Users/k-yoshimi/miniconda3/lib/libpython3.10.a \
  -DCMAKE_CXX_FLAGS="-isysroot $SDK -isystem $SDK/usr/include/c++/v1"
make _bse_solver && make install
cd ..
# Assertions on the installed tree:
LIBDIR=$(find "$PREFIX" -type d -name bse-python)
test -f "$LIBDIR/chiq/__init__.py"                     # chiq installed
test -f "$LIBDIR"/chiq/_bse_solver*.so                 # extension inside chiq
test -f "$LIBDIR/bse/__init__.py" && ! test -L "$LIBDIR/bse"   # real bse shim, NOT a symlink
test -f "$LIBDIR/bse_solver/__init__.py"               # top-level bse_solver shim
PYTHONPATH="$LIBDIR" python3 -c "from chiq import _bse_solver; import bse, bse_solver; print('legacy layout OK')"
echo "LEGACY INSTALL VERIFIED at $PREFIX"
```

- [ ] **Step 2: Run it to verify it fails**

Run: `bash tests/python/non-mpi/packaging/verify_legacy_install.sh`
Expected: FAIL — `CHIQ_WHEEL_BUILD` undefined is fine, but `python/CMakeLists.txt` still installs the `bse` **symlink** (so `! test -L "$LIBDIR/bse"` fails) and does not install the `bse`/`bse_solver` shim packages.

- [ ] **Step 3: Add the CMake option and dual layouts**

Top-level `CMakeLists.txt` — add the option near the top (after `project(ChiQ)`), and gate pybind11 acquisition + chiqvars:
```cmake
option(CHIQ_WHEEL_BUILD "Install layout for a pip/scikit-build-core wheel" OFF)
```
Change the pybind11 line:
```cmake
add_subdirectory(extern/pybind11)  # require cmake 3.4
```
to:
```cmake
if(CHIQ_WHEEL_BUILD)
  find_package(pybind11 CONFIG REQUIRED)   # sdists have no submodule; use the build-dep
else()
  add_subdirectory(extern/pybind11)
endif()
```
And skip the `chiqvars.sh` install in wheel mode — wrap the three chiqvars lines:
```cmake
if(NOT CHIQ_WHEEL_BUILD)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/chiqvars.sh.build.in ${CMAKE_CURRENT_BINARY_DIR}/chiqvars.sh @ONLY)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/chiqvars.sh.installed.in ${CMAKE_CURRENT_BINARY_DIR}/share/chiqvars.sh @ONLY)
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/share/chiqvars.sh DESTINATION share)
endif()
```

Rewrite `python/CMakeLists.txt` body (after the include guard) to:
```cmake
if(CHIQ_WHEEL_BUILD)
  # Wheel: pure-Python files come from scikit-build-core wheel.packages; the CLI
  # from pyproject entry points. Nothing to install here (the extension is
  # installed by src/CMakeLists.txt to DESTINATION chiq).
else()
  # Legacy layout.
  file(GLOB scripts LIST_DIRECTORIES false scripts/*.py)
  install(FILES ${scripts}
   PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ GROUP_EXECUTE GROUP_READ WORLD_EXECUTE WORLD_READ
   DESTINATION ${CMAKE_INSTALL_BINDIR})
  install(DIRECTORY package/chiq DESTINATION ${CMAKE_INSTALL_LIBDIR}/bse-python)
  # backward-compat shim packages (real dirs now, not a symlink)
  install(DIRECTORY package/bse        DESTINATION ${CMAKE_INSTALL_LIBDIR}/bse-python)
  install(DIRECTORY package/bse_solver DESTINATION ${CMAKE_INSTALL_LIBDIR}/bse-python)
  # remove any stale `bse` symlink from a previous install before the real dir lands
  install(CODE "file(REMOVE \"\$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}/bse-python/bse\")")
endif()
```
Note the `install(CODE file(REMOVE ...))` must run BEFORE the `install(DIRECTORY package/bse ...)` copies the real dir — CMake runs install rules in declaration order, so place the `file(REMOVE)` line ABOVE the `install(DIRECTORY package/bse ...)` line (reorder accordingly). Also remove the old `create_symlink` line entirely.

- [ ] **Step 4: Run the verification to confirm it passes**

Run: `bash tests/python/non-mpi/packaging/verify_legacy_install.sh`
Expected: prints `LEGACY INSTALL VERIFIED` — the installed `bse` is a real directory (not a symlink), `bse_solver` shim present, extension inside `chiq/`, all three importable.

- [ ] **Step 5: Commit**

```bash
git add CMakeLists.txt python/CMakeLists.txt tests/python/non-mpi/packaging/verify_legacy_install.sh
git commit -m "feat(build): CHIQ_WHEEL_BUILD option, pybind11 find_package, dual install layouts + shims"
```

---

### Task B3: `pyproject.toml` (scikit-build-core, deps, entry points, extras)

**Files:**
- Create: `pyproject.toml`
- Test: `tests/python/non-mpi/packaging/test_pyproject_metadata.py`

**Interfaces:**
- Consumes: `chiq.cli.<name>:main` (Plan A), `chiq.cli._deprecated.<name>_py` (Plan A), `CHIQ_WHEEL_BUILD` (B2).
- Produces: an installable/buildable project; `[project.scripts]` entry points.

- [ ] **Step 1: Write the failing test**

```python
# tests/python/non-mpi/packaging/test_pyproject_metadata.py
import os, sys
if sys.version_info >= (3, 11):
    import tomllib
    def _load(p):
        with open(p, "rb") as f: return tomllib.load(f)
else:
    import toml
    def _load(p): return toml.load(p)

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(__file__)))))

def test_pyproject_core_metadata():
    d = _load(os.path.join(ROOT, "pyproject.toml"))
    assert d["build-system"]["build-backend"] == "scikit_build_core.build"
    assert d["project"]["requires-python"] == ">=3.8"
    core = " ".join(d["project"]["dependencies"])
    for dep in ["numpy", "scipy", "more-itertools", "h5py", "toml"]:
        assert dep in core
    scripts = d["project"]["scripts"]
    for cmd in ["chiq_main", "chiq_post", "gen_qpath", "dcore_chiq"]:
        assert scripts[cmd] == f"chiq.cli.{cmd}:main"
    assert scripts["chiq_main.py"] == "chiq.cli._deprecated:chiq_main_py"
    extras = d["project"]["optional-dependencies"]
    assert "matplotlib" in " ".join(extras["plot"])
    assert "mpi4py" in " ".join(extras["mpi"])
    assert "dcore" in " ".join(extras["dcore"])
```

- [ ] **Step 2: Run test to verify it fails**

Run: `PYTHONPATH="$(pwd)/python/package" python3 -m pytest tests/python/non-mpi/packaging/test_pyproject_metadata.py -v`
Expected: FAIL — no `pyproject.toml`.

- [ ] **Step 3: Create `pyproject.toml`**

```toml
[build-system]
requires = ["scikit-build-core>=0.10,<1.0", "pybind11>=2.12,<3.0", "cmake>=3.15"]
build-backend = "scikit_build_core.build"

[project]
name = "chiq"
dynamic = ["version"]
description = "Solve the Bethe-Salpeter equation within DMFT; momentum-dependent susceptibilities."
readme = "README.md"
license = { file = "LICENSE" }
requires-python = ">=3.8"
dependencies = [
    "numpy>=1.23",
    "scipy",
    "more-itertools",
    "h5py",
    "toml",
]

[project.optional-dependencies]
plot = ["matplotlib"]
mpi = ["mpi4py"]
dcore = ["dcore"]
gpu = ["cupy"]
test = ["pytest"]

[project.scripts]
chiq_main = "chiq.cli.chiq_main:main"
chiq_post = "chiq.cli.chiq_post:main"
chiq_fft = "chiq.cli.chiq_fft:main"
gen_qpath = "chiq.cli.gen_qpath:main"
gen_allq = "chiq.cli.gen_allq:main"
calc_Iq = "chiq.cli.calc_Iq:main"
calc_Iq_scl = "chiq.cli.calc_Iq_scl:main"
plot_chiq_path = "chiq.cli.plot_chiq_path:main"
plot_Ir = "chiq.cli.plot_Ir:main"
eigenvec_viewer = "chiq.cli.eigenvec_viewer:main"
dcore_chiq = "chiq.cli.dcore_chiq:main"
"chiq_main.py" = "chiq.cli._deprecated:chiq_main_py"
"chiq_post.py" = "chiq.cli._deprecated:chiq_post_py"
"chiq_fft.py" = "chiq.cli._deprecated:chiq_fft_py"
"gen_qpath.py" = "chiq.cli._deprecated:gen_qpath_py"
"gen_allq.py" = "chiq.cli._deprecated:gen_allq_py"
"calc_Iq.py" = "chiq.cli._deprecated:calc_Iq_py"
"calc_Iq_scl.py" = "chiq.cli._deprecated:calc_Iq_scl_py"
"plot_chiq_path.py" = "chiq.cli._deprecated:plot_chiq_path_py"
"plot_Ir.py" = "chiq.cli._deprecated:plot_Ir_py"
"eigenvec_viewer.py" = "chiq.cli._deprecated:eigenvec_viewer_py"
"dcore_chiq.py" = "chiq.cli._deprecated:dcore_chiq_py"

[tool.scikit-build]
wheel.packages = ["python/package/chiq", "python/package/bse", "python/package/bse_solver"]
cmake.args = ["-DCHIQ_WHEEL_BUILD=ON", "-DTesting=OFF", "-DCMAKE_POLICY_VERSION_MINIMUM=3.5"]
sdist.include = ["src", "cmake", "extern", "tests", "python", "CMakeLists.txt", "LICENSE", "README.md"]

[tool.scikit-build.metadata.version]
provider = "scikit_build_core.metadata.regex"
input = "python/package/chiq/__init__.py"
regex = '__version__\s*=\s*["\'](?P<value>[^"\']+)["\']'
```

- [ ] **Step 4: Run test to verify it passes**

Run: `PYTHONPATH="$(pwd)/python/package" python3 -m pytest tests/python/non-mpi/packaging/test_pyproject_metadata.py -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add pyproject.toml tests/python/non-mpi/packaging/test_pyproject_metadata.py
git commit -m "feat(build): pyproject.toml (scikit-build-core, deps, entry points, extras)"
```

---

### Task B4: Clean-venv install tests (wheel, sdist, editable)

**Files:**
- Create: `tests/python/non-mpi/packaging/verify_pip_install.sh`

**Interfaces:**
- Consumes: `pyproject.toml` (B3), the CMake layouts (B1-B2).
- Produces: proof that `pip install .`, an sdist, and `pip install -e .` all work and import correctly in a clean venv.

- [ ] **Step 1: Write the failing verification script**

Create `tests/python/non-mpi/packaging/verify_pip_install.sh`:
```bash
#!/usr/bin/env bash
set -euo pipefail
REPO="$(cd "$(dirname "$0")/../../../.." && pwd)"
SDK=/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk
export CMAKE_ARGS="-DCMAKE_POLICY_VERSION_MINIMUM=3.5 \
  -DEIGEN3_INCLUDE_DIR=$(brew --prefix eigen)/include/eigen3 \
  -DPYTHON_INCLUDE_DIR=/Users/k-yoshimi/miniconda3/include/python3.10 \
  -DPYTHON_LIBRARY=/Users/k-yoshimi/miniconda3/lib/libpython3.10.a \
  -DCMAKE_CXX_FLAGS=-isysroot\\ $SDK\\ -isystem\\ $SDK/usr/include/c++/v1"

VENV=$(mktemp -d)/venv
python3 -m venv "$VENV"
# shellcheck disable=SC1091
source "$VENV/bin/activate"
pip install -U pip build

# --- wheel: build, install into the clean venv OUTSIDE the repo, import-check ---
WORK=$(mktemp -d); cp -R "$REPO" "$WORK/ChiQ"
cd "$WORK/ChiQ"
git submodule update --init >/dev/null 2>&1 || true   # legacy submodule (wheel uses find_package, but keep clean)
pip install .                                          # builds via scikit-build-core
cd /tmp                                                # import from OUTSIDE the source tree
python -c "import chiq, chiq._bse_solver; print('wheel import chiq OK')"
python -W ignore::FutureWarning -c "import bse, bse_solver; print('wheel shims OK')"
python -c "import sys; import chiq.g2scl_core; assert 'matplotlib' not in sys.modules; print('matplotlib optional OK')"
chiq_main --version                                   # console entry point on PATH
python - <<'PY'
import subprocess, sys
# deprecated .py alias prints a stderr notice
r = subprocess.run(["chiq_main.py", "--version"], capture_output=True, text=True)
assert "deprecated" in r.stderr.lower(), r.stderr
print("deprecated alias notice OK")
PY
# assert the installed wheel tree contains no test .h5 (artifact policy)
SITE=$(python -c "import chiq, os; print(os.path.dirname(os.path.dirname(chiq.__file__)))")
! find "$SITE/chiq" -name '*.h5' | grep -q . && echo "no test .h5 in wheel OK"
pip uninstall -y chiq

# --- sdist: build, inspect archive contents, install from it ---
cd "$WORK/ChiQ"
python -m build --sdist --outdir dist_sdist
SDIST=$(ls dist_sdist/*.tar.gz | head -1)
tar tzf "$SDIST" | grep -q 'src/bse_solver_pybind.cpp' && echo "sdist has C++ sources OK"
tar tzf "$SDIST" | grep -q 'tests/.*\.h5' && echo "sdist has test .h5 OK"
pip install "$SDIST"
cd /tmp && python -c "import chiq, chiq._bse_solver; print('sdist install OK')"
pip uninstall -y chiq

# --- editable ---
cd "$WORK/ChiQ"
pip install -e .
cd /tmp && python -c "import chiq, chiq._bse_solver; print('editable install OK')"
deactivate
echo "ALL PIP INSTALL PATHS VERIFIED"
```
Make it executable: `chmod +x tests/python/non-mpi/packaging/verify_pip_install.sh`.

- [ ] **Step 2: Run it to verify it fails / progresses**

Run: `bash tests/python/non-mpi/packaging/verify_pip_install.sh`
Expected initially: it may FAIL at `pip install .` if `CMAKE_ARGS` quoting or the scikit-build-core config needs adjustment. Debug iteratively: the most common issues are (a) the `-isysroot`/`-isystem` flags needing correct shell quoting inside `CMAKE_ARGS`; (b) `find_package(pybind11)` not finding the pip-installed pybind11 (scikit-build-core normally injects it — ensure `pybind11` is in `build-system.requires`); (c) `wheel.packages` paths. Fix `pyproject.toml`/the script until it prints `ALL PIP INSTALL PATHS VERIFIED`.

- [ ] **Step 3: (no separate implementation — the fixes land in `pyproject.toml`/CMake from B2-B3)**

Iterate on `pyproject.toml` / the CMake option until the script passes. If a real defect in B1-B3 is found, fix it there and note it.

- [ ] **Step 4: Run to verify it passes**

Run: `bash tests/python/non-mpi/packaging/verify_pip_install.sh`
Expected: prints `ALL PIP INSTALL PATHS VERIFIED`.

- [ ] **Step 5: Commit**

```bash
git add tests/python/non-mpi/packaging/verify_pip_install.sh pyproject.toml
git commit -m "test(build): clean-venv wheel/sdist/editable install verification"
```

---

### Task B5: CI job + documentation

**Files:**
- Modify: `.github/workflows/main.yml`
- Modify: `README.md`, `doc/install.rst`, `doc/tutorial/multiple_orb.rst`, `doc/tutorial/intersite_interactions.rst`, `doc/tutorial/single_orb.rst`
- Test: CI runs on push (the workflow itself is the test); docs are prose.

**Interfaces:**
- Consumes: everything above.
- Produces: a CI job exercising the pip path + docs describing it.

- [ ] **Step 1: Add a pip CI job to `.github/workflows/main.yml`**

Keep the existing `main` job (cmake→make→ctest→pytest). Add a sibling job (same `strategy.matrix.python-version: ['3.8','3.13']`) that installs system deps (`libeigen3-dev hdf5-tools openmpi-bin libopenmpi-dev`), then:
```yaml
      - name: pip build + install + test
        run: |
          python -m pip install -U pip
          pip install .[test,plot]
          python -c "import chiq, chiq._bse_solver; print('import OK')"
          python -W ignore::FutureWarning -c "import bse, bse_solver"
          chiq_main --version
          # run the non-mpi suite against the installed package (from outside the repo tree)
          cd /tmp && python -m pytest $GITHUB_WORKSPACE/tests/python/non-mpi -q -k "not mpi"
```
(On Ubuntu CI, `find_package(pybind11)` is satisfied by the `pybind11` build-requirement; no `CMAKE_ARGS` workaround is needed — those are local-macOS-only.) Add a lightweight 3.11 import-smoke job per the spec.

- [ ] **Step 2: Update the docs (spec §7-§8 inventory)**

- `README.md`: add a "Quick install (pip)" block (`pip install .`, `pip install -e .`, extras `pip install .[mpi,dcore]`) above the CMake section; change `chiq_main.py`/`chiq_post.py` usage examples to `chiq_main`/`chiq_post`; fix `sphinx-build -b html ../bse/doc html` → `../ChiQ/doc`.
- `doc/install.rst`: add the pip path; the **pip-vs-CMake decision matrix** (pip for laptops/venvs/CI; CMake for HPC/custom compilers); the **offline HPC** path (`pip install scikit-build-core pybind11 cmake` then `pip install . --no-build-isolation`); the cluster `mpi4py` note (`MPICC=mpicc pip install mpi4py --no-binary=mpi4py` after `module load`); the editable-install C++-rebuild caveat; the **unset-stale-`PYTHONPATH`** warning; change `chiq_main.py` → `chiq_main` (line ~74).
- `doc/tutorial/{multiple_orb,intersite_interactions,single_orb}.rst`: strip `.py` from every command example; fix the `chi_main.py` → `chiq_main` typo in `single_orb.rst`.

- [ ] **Step 3: Verify docs build**

Run: `sphinx-build -b html doc /tmp/built_doc 2>&1 | tail -5` (expect no new errors; needs `sphinx wild_sphinx_theme matplotlib`).

- [ ] **Step 4: Commit**

```bash
git add .github/workflows/main.yml README.md doc/install.rst doc/tutorial/*.rst
git commit -m "ci+docs: pip install job; document pip path, HPC offline, decision matrix; drop .py command names"
```

---

## Final verification

- [ ] `bash tests/python/non-mpi/packaging/verify_legacy_install.sh` → `LEGACY INSTALL VERIFIED`.
- [ ] `bash tests/python/non-mpi/packaging/verify_pip_install.sh` → `ALL PIP INSTALL PATHS VERIFIED`.
- [ ] Full non-mpi suite green against the source tree (with the `.so` copied into `python/package/chiq/` per B1 Step 5, or via an editable install).
- [ ] Legacy `cmake -DTesting=ON && make && source chiqvars.sh && ctest -V` still builds and the C++ tests pass (coexistence).
- [ ] `git grep -n "chiq_main.py" README.md doc/` returns only intentional references (deprecated-alias docs), not usage examples.

## Notes for the executor

- The macOS/clang-21 `CMAKE_ARGS` block (Global Constraints) is **local-only**; do not commit it into `pyproject.toml`. On Linux CI it is unnecessary.
- After Task B1 the built extension lives at `python/package/chiq/_bse_solver*.so` (gitignored). If the source-tree PYTHONPATH workflow misbehaves, prefer `pip install -e .` for local runs.
- If `pip install .` cannot build in this environment despite the `CMAKE_ARGS`, that is an environment/toolchain blocker (not a plan defect) — capture the exact error and escalate; the legacy `make install` path (Task B2 verification) still proves the packaging layout independently of pip.
