# Design: pip packaging via scikit-build-core (coexisting with the CMake workflow)

**Date:** 2026-07-12
**Status:** Ready for review (incorporates Codex + Antigravity design review rounds 1-3)
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
- `[tool.scikit-build]` config (exact):
  - **Pure-Python payload** — three package directories under `python/package/`:
    `wheel.packages = ["python/package/chiq", "python/package/bse", "python/package/bse_solver"]`.
    All three are packages (the `bse_solver` shim is a package, not a single `.py`, precisely
    so `wheel.packages` covers it with no single-module special-case). The §7 wheel test
    asserts the exact archive paths `chiq/__init__.py`, `bse/__init__.py`,
    `bse_solver/__init__.py`.
  - **`wheel.install-dir` is left UNSET** (default = wheel root). The extension is the only
    CMake-installed artifact in wheel mode and is installed with `install(TARGETS
    _bse_solver ... DESTINATION chiq)`, landing at `chiq/_bse_solver.<ext>` (no `chiq/chiq`
    duplication — the pure-Python `chiq/` comes from `wheel.packages`, the `.so` from this
    single CMake install rule). The §7 wheel test asserts `chiq/_bse_solver.<ext>` exists.
  - `cmake.args = ["-DCHIQ_WHEEL_BUILD=ON", "-DTesting=OFF", "-DCMAKE_POLICY_VERSION_MINIMUM=3.5"]`
    — a **new CMake option `CHIQ_WHEEL_BUILD`** (default `OFF`) selects wheel behavior: it
    installs *only* the `_bse_solver` target to `DESTINATION chiq` and **skips** the legacy
    scripts/`chiqvars.sh`/`bse-python`/symlink install rules (the wheel gets its CLI from
    entry points, §4; its pure-Python files from `wheel.packages`). When `OFF` (the legacy
    default) the current destinations are used, the extension installs to
    `lib/bse-python/chiq/`, and the two shims are installed into `lib/bse-python` (§3.3) so
    both layouts expose the same import surface. (The `CMAKE_POLICY_VERSION_MINIMUM=3.5` arg
    lets CMake 4.x configure the project, avoiding a build-dep upper bound on `cmake`.)
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
- **Top-level `bse_solver` shim — a package, not a single module**
  (`python/package/bse_solver/__init__.py`, `from chiq._bse_solver import *`) preserves
  `import bse_solver; bse_solver.BSESolver`. It is a *package* (directory + `__init__.py`)
  rather than a standalone `bse_solver.py` so it is covered by `wheel.packages` (which selects
  package directories, not single top-level modules) with no special single-file rule.
  Installed in **both** layouts (wheel: via `wheel.packages`; legacy CMake: into
  `lib/bse-python`, a new install rule — this closes the gap that the legacy path previously
  installed only `chiq`). Emits its own visible deprecation notice on first import (§10),
  independent of the `bse` shim's warning.
- **`bse` package — static forwarding shims, not a dynamic finder** (`python/package/bse/`):
  the previous "eager `sys.modules` registration" plan broke core installs (it would import
  `chiq.sumk_dft_chi`/`chiq.g2scl_core` and their optional deps), and a dynamic
  `MetaPathFinder` promising strict module-object identity is unsound (importlib can
  overwrite the canonical module's `__spec__`/`__path__`; reload/concurrency/cleanup are
  fragile) and breaks IDEs/linters/mypy. Instead, `bse/` is a **real package of thin static
  forwarding modules**, one per public submodule on an explicit **allowlist**
  (`h5bse`, `bse_toml`, `matrix_dict`, `point_group`, `sumk_dft_chi`, `g2scl_core`, `tools`,
  `index_pair`, `mpi`, and the `point_group_data` subpackage): each `bse/<mod>.py` does
  `from chiq.<mod> import *` (and `from chiq.<mod> import <mod> as ...` is not needed).
  - **Lazy by construction:** `import bse` runs only `bse/__init__.py` (which emits the
    deprecation notice and imports nothing heavy); `from bse import sumk_dft_chi` imports
    `bse/sumk_dft_chi.py` → `chiq.sumk_dft_chi` only then, so optional deps (`dcore`,
    `mpi4py`) are pulled only when that specific legacy submodule is used — never on a plain
    `import bse`.
  - **Identity contract (weakened, honest):** re-exported *classes/functions* are the same
    objects (`bse.h5bse.h5BSE is chiq.h5bse.h5BSE`), which is what legacy `from bse import
    h5bse; h5bse.h5BSE(...)` code needs; the *module objects* differ (`bse.h5bse is not
    chiq.h5bse`). This is deterministic, tooling-friendly, and needs no import hook.
  - `point_group_data` is mirrored as a real subpackage `bse/point_group_data/` with a thin
    forwarding module for **every** public module in `chiq/point_group_data/` — the full
    inventory is `__init__`, `base`, `C1`, `C2`, `D3`, `D4`, `D6`, `O`, `Oh` — so
    `import bse.point_group_data.<X>` works for all of them (tested, not just `C1`).
  - Supported/tested forms (§7): `import bse`, `import bse.h5bse`, `from bse import h5bse`,
    `import bse.point_group_data.C1`, mixed `chiq`/`bse` order.
  - Replaces the filesystem symlink (unreliable in wheels). **Legacy CMake install must first
    remove any existing `bse` symlink at the destination** *before* installing the real
    `bse/` package — otherwise `make install` on an upgrade could follow the old
    `bse -> chiq` symlink and overwrite `chiq/__init__.py`. Removal is done at **install
    time** via a portable `install(CODE "file(REMOVE ...)")` (not a shell `rm`, for
    Windows/portability), constructed to honor **`DESTDIR`** and `CMAKE_INSTALL_PREFIX`
    (staged/packaged installs) and to target *only* `<dest>/lib/bse-python/bse`. Both layouts
    then ship the same real shim and one deprecation behavior.
  - **Static-shim caveat (documented):** `from chiq.<mod> import *` forwards only the public
    names (respecting each module's `__all__`; underscore-prefixed names are not re-exported).
    The allowlisted modules' public surfaces are the compatibility contract; a legacy caller
    reaching a private `bse.<mod>._x` name is out of scope (none are known to).

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
- **Entry points** under `[project.scripts]`: `chiq_main = "chiq.cli.chiq_main:main"`, etc. —
  one per command. Optional-dependency commands (`dcore_chiq` needs `dcore`; `plot_chiq_path`,
  `plot_Ir` need `matplotlib`) still install as entry points, but each `main()` **parses
  `--help`/`--version` (argparse) BEFORE importing its optional dep** — the heavy import
  happens only after arg parsing, inside the code path that does real work. So
  `dcore_chiq --help` / `--version` succeed on a **core** install (tested in the core-only
  env, §7), while actually running it without the extra yields a clear
  `"<command> requires the '<extra>' extra: pip install chiq[<extra>]"` error.
- **Legacy compatibility.** `python/scripts/<name>.py` are kept as **thin wrappers**
  (`from chiq.cli.<name> import main; main()`) so the direct-CMake path installs them to
  `bin` exactly as today. Additionally, deprecated `<name>.py` console-script aliases are
  shipped for one release as **quoted** script keys (dotted names are invalid unquoted TOML
  keys), e.g. `"chiq_main.py" = "chiq.cli._deprecated:chiq_main_py"`; each prints a
  deprecation line to `stderr` before delegating (so shell users see it, §10). A §7 test
  asserts the generated executable has the exact legacy name.

## 5. Optional-dependency handling (module-scope import audit)

Confirmed by grep of the actual source:

- `chiq.sumk_dft_chi`, `chiq.dcore_chiq_worker`, `chiq.cli.dcore_chiq` import **dcore**
  (and `mpi4py`) at module scope → gated behind the `dcore`/`mpi` extras; these are only
  reached by the χ₀ preprocessing command, imported lazily in that command's `main()`.
- `chiq.g2scl_core` imports **matplotlib** at module scope (`matplotlib.use('Agg')`). To keep
  matplotlib **optional** (headless HPC nodes should not need it just to solve), this design
  makes that import **lazy**: move the `matplotlib`/`pyplot`/`cm` imports out of
  `chiq.g2scl_core` module scope into the plotting functions that use them (a small, behavior-
  preserving refactor), and likewise keep the plotting CLI modules (`plot_chiq_path`,
  `plot_Ir`) importing matplotlib lazily inside `main()`. matplotlib then stays a `plot`
  extra; importing `chiq.g2scl_core` (or running a non-plot command) no longer requires it,
  and invoking a plot command without `[plot]` yields a clear missing-extra error.
- `chiq._mpi` imports `mpi4py`, but `chiq.mpi` falls back to `chiq._no_mpi` when absent, so
  core commands (`chiq_main` etc.) run without `mpi4py`.

## 6. Dependencies and Python matrix (pinned)

- **Python:** `requires-python = ">=3.8"` (kept for ohtaka). CI runs the full
  build+test on **3.8 and 3.13**; a note states only these endpoints are continuously
  validated (a 3.11 wheel/import smoke job is added to catch mid-range breakage).
- **Core runtime deps:** `numpy>=1.23`, `scipy`, `more-itertools`, `h5py`, `toml`.
  (`toml` is kept because `bse_toml.py` both reads *and* writes TOML — stdlib `tomllib`/
  `tomli` are read-only, so a straight swap is not possible; migrating the writer is a
  separate future cleanup, out of scope here.)
- **Optional extras:** `plot` → `matplotlib` (lazy, §5); `mpi` → `mpi4py`; `dcore` → `dcore`
  (χ₀ preprocessing); `gpu` → `cupy` (future); `test` → `pytest`.
- **Build deps:** `scikit-build-core>=0.10,<1.0`, `pybind11>=2.12,<3.0`, `cmake>=3.15`
  (no `cmake` upper bound — CMake 4.x is supported via the
  `-DCMAKE_POLICY_VERSION_MINIMUM=3.5` arg in §3.1). The **sdist build is exercised under
  Python 3.8 in CI** (§7) so build-dep drift that drops 3.8 is caught early rather than only
  by `requires-python` advertising.

## 7. Testing

- **Wheel/sdist install tests (the crux of "pip works"):** build both an sdist and a wheel;
  `pip install` each into a **clean venv outside the repo** and smoke-test: `import chiq`,
  `import chiq._bse_solver` (extension loads), the `bse`/`bse_solver` shims incl. nested and
  mixed-order forms (§3.3), `point_group_data` load, and each **core** console command's
  `--version`/`--help` (`chiq_main`, `chiq_post`, `chiq_fft`, `gen_qpath`, `gen_allq`,
  `calc_Iq`, `calc_Iq_scl`, `plot_chiq_path`, `plot_Ir`, `eigenvec_viewer`). The
  extra-gated commands' **`--help`/`--version` are tested in the CORE-only env** (they must
  succeed there because help/version are parsed before any optional import, §4); their
  **real execution** without the extra is tested to emit the actionable missing-extra message
  (via a `try/except ImportError` in `main()`, not a raw traceback), and a separate `[dcore]`
  env exercises `dcore_chiq` end-to-end. Run under 3.8 and 3.13.
- **sdist content assertion:** unpack the **sdist archive itself** (not the installed env,
  which builds a wheel) and assert it contains everything needed to configure+compile in a
  clean checkout-less environment — the required CMake files (`CMakeLists.txt`,
  `src/CMakeLists.txt`, `cmake/`), C++ sources/headers, all Python sources
  (`chiq/`, `bse/`, `bse_solver.py`, `chiq/cli/`), `point_group_data`, test `.h5` fixtures,
  and `LICENSE`/`README.md`. Additionally build the sdist in an environment **without the
  git submodule** to confirm the `find_package(pybind11)` path (§3.2) works. Assert the
  built **wheel** archive contains `chiq/` (with `chiq/_bse_solver.<ext>` + `point_group_data`),
  `bse/`, `bse_solver.py`, and **no** test `.h5`, `lib/bse-python`, raw scripts, or
  `chiqvars.sh`.
- **Editable install:** `pip install -e .` then `import chiq._bse_solver` from a directory
  **outside** the repo (so in-tree `.so`/source can't mask a broken editable install);
  document that C++ source edits require an explicit `pip install -e . --no-build-isolation`
  rebuild (scikit-build-core does not auto-recompile on import).
- **Direct-CMake path unchanged:** a CI job still runs `cmake -DTesting=ON && make &&
  source chiqvars.sh && ctest -V && pytest`, proving the legacy workflow still works and
  installs the extension at the same `chiq._bse_solver` location. Run the legacy import
  checks from a directory **outside the repo** and assert each imported module's `__file__`
  lies under the temporary install prefix (not the source tree), so the test validates the
  *installed* shims/extension rather than in-tree files. Also verify the old `bse` symlink is
  gone and `bse/` is the real shim package (the §3.3 symlink-removal rule worked).
- **Shim tests:** `bse`/`bse_solver` import-form/identity/import-order tests. `bse` and
  `bse_solver` are **independent** deprecated surfaces — each emits its own notice once per
  process (a test asserts `import bse` warns once and `import bse_solver` warns once,
  independently).

## 8. Files touched / added

- Add: `pyproject.toml`; `python/package/bse_solver/__init__.py` (top-level shim package);
  `python/package/bse/` (static forwarding-shim package: `__init__.py` + one thin
  `from chiq.<mod> import *` module per allowlisted submodule + `point_group_data/` mirror);
  `python/package/chiq/cli/` (one module per command + `_deprecated.py` wrappers).
- Edit: `src/bse_solver_pybind.cpp` (module rename); top-level `CMakeLists.txt` +
  `src/CMakeLists.txt` + `python/CMakeLists.txt` (add `CHIQ_WHEEL_BUILD`, pybind11 selection,
  install destinations/shims for both layouts); `python/package/chiq/solver/cpp.py`
  (`from chiq import _bse_solver`).
- Edit: `python/scripts/*.py` → thin wrappers delegating to `chiq.cli.*`.
- Edit: `.github/workflows/main.yml` (pip job + install/sdist/editable tests + 3.11 smoke).
- Edit docs: `README.md` — add a "Quick install (pip)" subsection (`pip install .`) above the
  CMake instructions; update the "How to use" examples from `chiq_main.py`/`chiq_post.py` to
  the entry-point names `chiq_main`/`chiq_post`; fix the stale doc-build path
  `sphinx-build -b html ../bse/doc html` → `../ChiQ/doc`. `doc/install.rst` — pip vs CMake
  decision matrix (§9), offline `--no-build-isolation` path (§9), cluster `mpi4py` build
  (§9), editable-install C++ rebuild caveat, and a **warning to unset a stale `PYTHONPATH`**
  from a prior `chiqvars.sh` when switching to the pip install (else the old tree masks the
  pip package / a stale `bse_solver` shadows the shim).

## 9. HPC / offline build (ohtaka)

- Default `pip install .` uses PEP 517 **build isolation**, which downloads the build
  requirements (`scikit-build-core`, `pybind11`, `cmake`) from PyPI — this **fails on
  air-gapped cluster nodes**. `doc/install.rst` documents the offline path:
  `pip install scikit-build-core pybind11 cmake` (or `module load` them) into the target
  env, then `pip install . --no-build-isolation`.
- Compiler / Eigen selection: cluster users pass toolchain and Eigen hints via
  `CMAKE_ARGS="-DEIGEN3_DIR=... -DCMAKE_CXX_COMPILER=..."` env var (scikit-build-core forwards
  it). Documented alongside.
- **`mpi` extra on clusters:** `pip install chiq[mpi]` builds `mpi4py`, which must link the
  cluster's MPI (not a generic host MPI). `doc/install.rst` instructs users to `module load`
  their MPI and pre-install with the cluster wrapper, e.g.
  `MPICC=mpicc pip install mpi4py --no-binary=mpi4py`, *before* installing chiq.
- **pip-vs-CMake decision matrix** (in `doc/install.rst`): recommend **pip** for standard
  users / laptops / virtualenvs / CI; recommend the **CMake + chiqvars.sh** workflow for HPC
  admins or anyone needing custom compilers, custom Eigen search paths, or specialized build
  flags. This is why the CMake path is retained, not merely tolerated.
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

## 11. Review resolution (Codex + Antigravity)

### Round 3

- **`bse_solver.py` single-module wheel inclusion undefined (Codex must-fix):** the top-level
  shim is now a **package** `bse_solver/__init__.py`, covered by `wheel.packages` like the
  others — no single-file mechanism needed (§3.1, §3.3, §8).
- **§8 still said "MetaPathFinder" (Codex should-fix):** file inventory corrected to the static
  forwarding-shim package (§8).
- **§4/§7 `dcore_chiq --help` env conflict (Codex should-fix):** `--help`/`--version` tested in
  the **core-only** env; missing-extra message and end-to-end tested with `[dcore]` (§7).
- **Symlink removal portability/DESTDIR (Codex + Antigravity):** removal via
  `install(CODE "file(REMOVE ...)")` honoring `DESTDIR`/prefix, targeting only
  `lib/bse-python/bse` (§3.3).
- **`point_group_data` nested inventory (Codex should-fix):** full module list enumerated
  (`base,C1,C2,D3,D4,D6,O,Oh`) and all tested (§3.3).
- **Per-shim warning semantics (Codex should-fix):** `bse` and `bse_solver` warn independently,
  once each per process (§7).
- **`import *` only forwards public names (Codex risk):** documented static-shim caveat (§3.3).
- **README stale (`chiq_main.py`, `../bse/doc`) + no pip mention (Antigravity must-fix-docs):**
  captured as required doc edits in §8 (README pip quick-install, `chiq_main` command names,
  fixed doc-build path) + PYTHONPATH-unset warning — applied during implementation, not now
  (design HARD-GATE).
- Unpinned runtime deps vs 3.8 (Codex risk): accepted — endpoint CI catches drift; hard pins
  would over-constrain research users. Documented as a monitored risk.

### Round 2

- **`bse` MetaPathFinder unsound + breaks tooling (Codex must-fix, Antigravity risk):**
  replaced with a **static forwarding-shim package** (allowlist of thin `bse/<mod>.py`
  re-export modules); lazy per-submodule, deterministic, IDE/mypy-friendly; identity contract
  weakened to class/function identity (honest) rather than module-object identity (§3.3).
- **Wheel shim inclusion + `wheel.install-dir` ambiguity (Codex must-fix):** `wheel.packages`
  lists `chiq` and `bse`; `bse_solver.py` included as a top-level module; `wheel.install-dir`
  left unset with the extension installed to `DESTINATION chiq`; exact archive paths asserted
  (§3.1, §7).
- **CMake symlink-follow overwrites `chiq/__init__.py` on upgrade (Antigravity must-fix):**
  legacy install `rm -f`s the old `bse` symlink before installing the real `bse/` package
  (§3.3).
- **matplotlib bloats headless nodes (Antigravity should-fix):** made lazy in `chiq.g2scl_core`
  + plot CLI, kept as the `plot` extra rather than promoted to core (§5–6).
- **`cmake` unbounded / 3.8 (Codex should-fix):** no `cmake` upper bound; CMake 4.x handled via
  `-DCMAKE_POLICY_VERSION_MINIMUM=3.5`; sdist build tested on 3.8 (§3.1, §6).
- **Optional-CLI `--help` on core (Codex should-fix):** `--help`/`--version` parsed before the
  optional import; `dcore_chiq --help` tested in the core env (§4, §7).
- **Dotted script keys / legacy-origin / sdist allowlist (Codex should-fix):** quoted
  `"chiq_main.py"` keys; legacy import tests assert install-prefix origin; sdist asserted to
  build submodule-free with a full file allowlist (§4, §7).
- **mpi4py cluster linking + pip-vs-cmake guidance (Antigravity should-fix):** documented
  `MPICC=... pip install mpi4py --no-binary` and a decision matrix in `doc/install.rst` (§9).
- Editable-install rebuild caveat (both, open question): documented — scikit-build-core does
  not auto-recompile on import; C++ edits need an explicit `-e .` rebuild (§7).

### Round 1

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
