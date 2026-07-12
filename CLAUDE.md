# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What ChiQ is

ChiQ solves the Bethe-Salpeter equation (BSE) within DMFT and computes momentum-dependent susceptibilities χ(q, iν_m) for correlated electron systems. It runs as a post-processing step of [DCore](https://github.com/issp-center-dev/DCore) and supports several approximation schemes: BSE, SCL, RPA, and RRPA (defined in `doc/algorithms.rst`).

## Build

CMake-based, C++17 + Eigen3 (modern Eigen 3.4+ requires ≥ C++14), with a pybind11 submodule (`extern/pybind11` — run `git submodule init && git submodule update` after cloning). Build out of source:

```bash
mkdir build && cd build
cmake .. -DTesting=ON -DCMAKE_INSTALL_PREFIX=$HOME/local
make
make install   # optional
```

- If Eigen3 isn't found: `-DEIGEN3_DIR=/path/to/eigen3-prefix`.
- After building, `source ./chiqvars.sh` (generated in the build dir) to set PATH/PYTHONPATH for running scripts and tests against the build tree.

## Tests

- **C++ tests** (googletest, bundled in `tests/c++/`): run `ctest -V` (or `make test`) in the build directory. Requires `-DTesting=ON`.
- **Python tests**: run `pytest` from the build directory (CMake copies `tests/python/` there and writes a `pytest.ini` that puts `python/package` on `pythonpath`). Requires the built `bse_solver` module to be importable — source `chiqvars.sh` first.
  - Single test: `pytest tests/python/non-mpi/bsetool_BSE/test_bsetool_BSE.py`
  - MPI tests live in `tests/python/mpi/`, non-MPI in `tests/python/non-mpi/`. Test dirs are named per scheme (`bsetool_BSE`, `bsetool_RPA`, `bsetool_SCL`, `bsetool_RRPA`, ...), each with its own `bse.toml` input and reference data.
- CI (`.github/workflows/main.yml`) mirrors this exactly: cmake → make → `source ./chiqvars.sh` → `ctest -V` → `pytest`. Python deps used in CI: `mpi4py dcore more-itertools pytest`; system deps: `hdf5-tools libeigen3-dev openmpi`.

## Documentation

Sphinx docs in `doc/` (needs `sphinx wild_sphinx_theme matplotlib`):
```bash
sphinx-build -b html doc built_doc
```
Docs deploy automatically to gh-pages on pushes to `main`/`develop` and `*-autodoc` branches.

## Architecture

Two layers connected through pybind11:

1. **C++ core** (`src/`): the numerical BSE solver. `bse.hpp` holds the solver logic, `block_matrix.hpp` a block-diagonal matrix type. Only a Python module is built (`bse_solver` via `bse_solver_pybind.cpp`) — the standalone C++ executable is commented out in `src/CMakeLists.txt`.

2. **Python layer** (`python/`):
   - `python/package/chiq/` — the `chiq` package: HDF5 I/O (`h5bse.py`, all data flows through an HDF5 file produced by DCore), TOML input parsing (`bse_toml.py`), lattice Green's function / χ₀ machinery (`sumk_dft_chi.py`), SCL core (`g2scl_core.py`), point-group symmetry data (`point_group.py`, `point_group_data/`), and an MPI shim (`mpi.py` falls back to `_no_mpi.py` when mpi4py is absent).
   - `python/package/bse` is a **symlink** to `chiq` kept for backward compatibility (the install step also creates it) — never edit files through the `bse` path.
   - `python/package/chiq/solver/` — the solver backend switch. `get_solver(backend, ...)` returns a `set`/`calc`/`get` facade for `backend ∈ {"cpp", "numpy", "cupy"}` (selected via `[chiq_main] backend`, default `"cpp"`). `cpp.py` wraps the pybind `bse_solver` module; `numpy.py`+`kernels.py`+`layout.py` are a pure-NumPy reimplementation that reproduces the C++ results (validated by `tests/python/non-mpi/solver/test_backend_agreement.py`); `cupy.py` is a skeleton for a future GPU backend.
   - `python/scripts/` — user-facing CLI entry points, installed to `bin`. Main pipeline: `gen_qpath.py` (q-path generation) → `dcore_chiq.py` (compute χ₀ and vertex from DCore output) → `chiq_main.py` (MPI-parallel driver; constructs the solver via `chiq.solver.get_solver`) → `chiq_post.py` (eigenvalue/matrix-element post-processing; imports `chiq_main` from the same directory) → `plot_chiq_path.py`.

The end-to-end DCore → ChiQ workflow is demonstrated in `examples/*/run.sh` (e.g. `examples/square_bse/`), and input is a TOML file with `[chiq_common]`, `[chiq_main]`, `[chiq_post]` sections.

## Conventions

- Version string lives in `python/package/chiq/__init__.py` (`__version__`).
- CI tests Python 3.8 and 3.13 — keep the Python code compatible with 3.8 (no 3.9+-only syntax).
