# Pip Packaging Verification Remediation Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task by task.

**Goal:** Make ChiQ's wheel, sdist, editable, and legacy packaging checks prove artifact isolation, archive safety, CLI/shim compatibility, and the NumPy/pybind11 ABI contract on Python 3.8 and 3.13.

**Architecture:** A small Python verification library builds a descriptor-safe tracked-file snapshot, validates sdists and wheels as hostile archives, and performs origin/CLI/ABI runtime checks. Shell entry points orchestrate fresh virtual environments without reusing artifacts; CMake owns the single native extension and narrowly removes obsolete legacy binaries. CI invokes each authoritative verifier once per endpoint and preserves diagnostic manifests and logs.

**Tech Stack:** Python 3.8-compatible stdlib, pytest, packaging, Bash, CMake, scikit-build-core, pybind11, NumPy, GitHub Actions.

---

## Ground rules

- Work on `design/pip-packaging`; do not switch branches or create another worktree for this continuation.
- Preserve all existing untracked scientific-test outputs. Never add them to a commit.
- Use the local Apple-clang flags from `AGENTS.md` before every local CMake or pip build. Never put those machine paths in `pyproject.toml`.
- Add no Python 3.9+ syntax. In particular, use `typing.List`/`typing.Optional`, not built-in generic aliases or `X | None`.
- Run each red test before its implementation and record the expected failure named below.
- Commit only the files listed by each task. Run `git diff --cached --check` before every commit.
- The addendum supersedes the base design's DCore numerical end-to-end wording: DCore acceptance is the import boundary plus a nonexistent-input invocation, not a scientific run.

## Shared acceptance constants

Create one manifest rather than duplicating names across tests and scripts:

```python
CANONICAL_COMMANDS = (
    "chiq_main", "chiq_post", "chiq_fft", "gen_qpath", "gen_allq",
    "calc_Iq", "calc_Iq_scl", "plot_chiq_path", "plot_Ir",
    "eigenvec_viewer", "dcore_chiq",
)
DEPRECATED_COMMANDS = tuple(name + ".py" for name in CANONICAL_COMMANDS)
SHIM_MODULES = (
    "bse_toml", "g2scl_core", "h5bse", "index_pair", "matrix_dict",
    "mpi", "point_group", "sumk_dft_chi", "tools",
)
POINT_GROUP_MODULES = ("C1", "C2", "D3", "D4", "D6", "O", "Oh", "base")
NATIVE_SUFFIXES = (".so", ".pyd", ".dylib")
```

The implementation may place these in `tests/python/non-mpi/packaging/contracts.py`; production code must not import from tests.

### Task 1: Lock the metadata, CLI, and compatibility contracts

**Files:**

- Create: `tests/python/non-mpi/packaging/contracts.py`
- Create: `python/package/chiq/cli/_common.py`
- Modify: `python/package/chiq/cli/*.py` for the 11 canonical parsers
- Modify: `python/package/chiq/cli/_deprecated.py`
- Modify: `python/package/bse/__init__.py`
- Modify: `python/package/bse_solver/__init__.py`
- Modify: `tests/python/non-mpi/packaging/test_cli_package.py`
- Modify: `tests/python/non-mpi/packaging/test_bse_shim.py`
- Modify: `tests/python/non-mpi/packaging/test_bse_solver_shim.py`
- Modify: `tests/python/non-mpi/packaging/test_pyproject_metadata.py`

- [ ] **Step 1: Write failing contract tests.**

  Add parameterized subprocess tests that run every canonical command's `--help` and `--version`, assert exact stdout `ChiQ version <chiq.__version__>`, and assert no traceback. Run with a scrubbed `PYTHONPATH` from a directory outside the checkout. Keep valid existing parser invocations in the existing integration suite.

  Add alias tests that assert every `.py` entry point writes a deprecation notice directly to stderr and forwards the original arguments. Add fresh-process warning tests with `warnings.catch_warnings(record=True)`:

  ```python
  import importlib
  import warnings

  with warnings.catch_warnings(record=True) as caught:
      warnings.simplefilter("always")
      importlib.import_module("bse")
      importlib.import_module("bse.tools")
      importlib.import_module("bse.point_group.C1")
  assert [w.category for w in caught].count(FutureWarning) == 1
  ```

  Separately prove that importing `bse_solver` after `bse` emits its own single `FutureWarning`, that forwarding module objects are distinct, and that named public objects are identical. Do not compare module identity.

  Extend metadata tests to require exactly 22 scripts, retain the three `wheel.packages`, retain dynamic version sourcing from `chiq.__version__`, and reject any local macOS paths.

- [ ] **Step 2: Run the focused tests and confirm red.**

  ```bash
  PYTHONPATH="$PWD/python/package:$PWD/build/src" python3 -m pytest \
    tests/python/non-mpi/packaging/test_cli_package.py \
    tests/python/non-mpi/packaging/test_bse_shim.py \
    tests/python/non-mpi/packaging/test_bse_solver_shim.py \
    tests/python/non-mpi/packaging/test_pyproject_metadata.py -q
  ```

  Expected failures: canonical parsers other than the two existing exceptions reject `--version`; warning-count and complete manifest assertions are absent or fail.

- [ ] **Step 3: Add one production version helper and use it everywhere.**

  Implement `python/package/chiq/cli/_common.py` exactly around the public version source:

  ```python
  from chiq import __version__

  VERSION_MESSAGE = "ChiQ version %s" % __version__

  def add_version_argument(parser):
      parser.add_argument("--version", action="version", version=VERSION_MESSAGE)
  ```

  Import and call `add_version_argument(parser)` after constructing each parser. Remove duplicate version declarations. Do not hard-code `prog="dcore_chiq.py"`; argparse must derive the invoked name from `argv[0]`.

  Keep plot imports after argument parsing. Only `plot_chiq_path` and `plot_Ir` need a valid minimal-fixture test that reaches the missing-matplotlib gate and prints `pip install chiq[plot]` without a traceback. `eigenvec_viewer` receives ordinary core-only help/version coverage.

- [ ] **Step 4: Make warning behavior explicit.**

  Compatibility packages use standard warnings:

  ```python
  warnings.warn(message, FutureWarning, stacklevel=2)
  ```

  Ensure `bse` owns one process-global warning guard shared by its nested modules, while `bse_solver` owns an independent guard. Keep `_deprecated.py` as direct stderr output because command aliases are executable wrappers, not import shims.

- [ ] **Step 5: Re-run focused and existing non-MPI tests.**

  ```bash
  PYTHONPATH="$PWD/python/package:$PWD/build/src" PATH="$PWD/python/scripts:$PATH" \
    python3 -m pytest tests/python/non-mpi/packaging -q
  cd build && PYTHONPATH="$OLDPWD/python/package:$OLDPWD/build/src" \
    PATH="$OLDPWD/python/scripts:$PATH" python3 -m pytest ../tests/python/non-mpi -q
  ```

  Expected: all tests pass; no parser defaults or numerical outputs change.

- [ ] **Step 6: Commit.**

  ```bash
  git add python/package/chiq/cli python/package/bse/__init__.py \
    python/package/bse_solver/__init__.py tests/python/non-mpi/packaging/contracts.py \
    tests/python/non-mpi/packaging/test_cli_package.py \
    tests/python/non-mpi/packaging/test_bse_shim.py \
    tests/python/non-mpi/packaging/test_bse_solver_shim.py \
    tests/python/non-mpi/packaging/test_pyproject_metadata.py
  git diff --cached --check
  git commit -m "test(packaging): lock CLI and shim contracts"
  ```

### Task 2: Build a race-resistant tracked-source snapshot

**Files:**

- Create: `tests/python/non-mpi/packaging/source_snapshot.py`
- Create: `tests/python/non-mpi/packaging/test_source_snapshot.py`

- [ ] **Step 1: Write adversarial tests first.**

  Build temporary Git repositories and cover: ordinary `100644`/`100755` files, unstaged tracked bytes, staged new files, excluded untracked warnings, `160000` gitlink omission, unmerged entries, tracked symlinks, non-regular leaves, executable-mode mismatch, symlink ancestors, leaf replacement between `lstat` and `open`, ancestor replacement during traversal, unsupported `O_NOFOLLOW`, and prohibited ignored cache/native files under `python/package`.

  Monkeypatch the file-open seam to replace an ancestor after its descriptor is opened; assert the copy either uses the already-open safe descriptor or fails closed. Add normalized NFC + `casefold()` collision tests even on a case-sensitive filesystem.

- [ ] **Step 2: Run the tests and confirm red.**

  ```bash
  python3 -m pytest tests/python/non-mpi/packaging/test_source_snapshot.py -q
  ```

  Expected failure: `ModuleNotFoundError: No module named 'source_snapshot'`.

- [ ] **Step 3: Implement index parsing and descriptor-relative copying.**

  Parse `git ls-files --stage -z` without shell splitting. Reject stages other than zero and modes outside `100644`, `100755`, and skipped `160000`. Reject all `120000` entries. Walk each path from an open checkout-root descriptor using `os.open(component, O_DIRECTORY | O_NOFOLLOW, dir_fd=parent_fd)` for ancestors and `os.open(leaf, O_RDONLY | O_NOFOLLOW, dir_fd=parent_fd)` for leaves. Compare `os.fstat()` to descriptor-relative `os.stat(..., follow_symlinks=False)` and index mode before copying bytes from the open leaf.

  Required public seam:

  ```python
  def create_snapshot(repo, destination, stderr):
      """Copy validated tracked worktree bytes and return the manifest dict."""
  ```

  Fail closed if `dir_fd`, `O_NOFOLLOW`, or no-follow stat support is unavailable. Create the destination with mode `0700`; open output files exclusively and set only `0644`/`0755`. Return `HEAD`, `git diff HEAD --stat`, and sorted `{path, index_mode, sha256}` rows to the caller, which writes `snapshot-manifest.json` beside the retained artifacts—not inside the build snapshot.

  Warn for non-ignored untracked files from `git status --porcelain=v1 -z --untracked-files=all`. Independently scan ignored and untracked paths under package roots and fail for `__pycache__`, `*.pyc`, `*.so`, `*.pyd`, or `*.dylib`. This makes a tracked native artifact a hard failure and prevents it from being mistaken for the CMake output.

- [ ] **Step 4: Re-run tests and inspect the real manifest.**

  ```bash
  python3 -m pytest tests/python/non-mpi/packaging/test_source_snapshot.py -q
  tmp=$(mktemp -d)
  python3 tests/python/non-mpi/packaging/source_snapshot.py \
    "$PWD" "$tmp/source" --manifest "$tmp/snapshot-manifest.json"
  python3 -m json.tool "$tmp/snapshot-manifest.json"
  test ! -e "$tmp/source/.git"
  test ! -e "$tmp/source/extern/pybind11/CMakeLists.txt"
  ```

  Expected: tests pass; the real manifest omits the gitlink contents and reports the existing excluded untracked scientific outputs without copying them.

- [ ] **Step 5: Commit.**

  ```bash
  git add tests/python/non-mpi/packaging/source_snapshot.py \
    tests/python/non-mpi/packaging/test_source_snapshot.py
  git diff --cached --check
  git commit -m "test(packaging): add safe tracked-source snapshots"
  ```

### Task 3: Validate and safely extract sdists

**Files:**

- Create: `tests/python/non-mpi/packaging/artifact_inspector.py`
- Create: `tests/python/non-mpi/packaging/test_sdist_inspector.py`

- [ ] **Step 1: Create the hostile tar corpus.**

  Generate archives in memory for traversal, absolute/drive/UNC/backslash paths, NULs, repeated separators, `.`/`..`, Unicode/casefold aliases, duplicate entries, file/directory and ancestor-file conflicts, symlink/hardlink/device/FIFO/sparse members, >50,000 members, >256 MiB logical files, >2 GiB logical totals, >200 aggregate ratio, and archives >512 MiB. Cover an omitted explicit root directory, owner-executable permission preservation, non-owner mode stripping, and destination containment when `/tmp` resolves through `/private/tmp` on macOS.

- [ ] **Step 2: Run the tests and confirm red.**

  ```bash
  python3 -m pytest tests/python/non-mpi/packaging/test_sdist_inspector.py -q
  ```

  Expected failure: inspector imports or validation entry points are missing.

- [ ] **Step 3: Implement graph validation before writes.**

  Define hard limits and one canonical path function:

  ```python
  MAX_ARCHIVE = 512 * 1024 * 1024
  MAX_MEMBERS = 50000
  MAX_FILE = 256 * 1024 * 1024
  MAX_TOTAL = 2 * 1024 * 1024 * 1024
  MAX_RATIO = 200

  def validate_sdist(path):
      """Return a validated manifest without extracting any member."""

  def extract_sdist(path, destination):
      """Validate the complete graph, then manually extract regular files."""
  ```

  Normalize with `PurePosixPath`, NFC, and casefold collision keys. Require one normalized `<name>-<version>` top-level component and compare it later with metadata. Accept only regular files and directories; explicitly reject tar sparse indicators. Resolve both destination root and every parent before containment checks. Never use `extractall`; create files with exclusive mode, copy bounded streams, then apply `0755` only when the source owner-execute bit was set and `0644` otherwise. Directories become `0755`; never restore owner or timestamps.

- [ ] **Step 4: Add sdist content and metadata agreement checks.**

  Require `pyproject.toml`, `LICENSE`, `README.md`, root/python/src CMake files, C++ sources and headers, CMake helpers, all three package trees, CLI and point-group modules, and test HDF5 fixtures. Reject caches, installed layouts, build trees, and all native suffixes.

  Parse `PKG-INFO` after safe extraction and require normalized project name/version agreement among archive filename, top-level directory, `PKG-INFO`, and the dynamic version provider declared in `pyproject.toml` (`chiq.__version__`).

- [ ] **Step 5: Run tests and commit.**

  ```bash
  python3 -m pytest tests/python/non-mpi/packaging/test_sdist_inspector.py -q
  git add tests/python/non-mpi/packaging/artifact_inspector.py \
    tests/python/non-mpi/packaging/test_sdist_inspector.py
  git diff --cached --check
  git commit -m "test(packaging): validate and extract sdists safely"
  ```

### Task 4: Validate wheel ZIP, metadata, and RECORD contracts

**Files:**

- Modify: `tests/python/non-mpi/packaging/artifact_inspector.py`
- Create: `tests/python/non-mpi/packaging/test_wheel_inspector.py`

- [ ] **Step 1: Write malformed-wheel tests.**

  Cover encrypted/multi-disk/overlapping/out-of-range entries, bad CRC, path aliases and collisions, per-entry and aggregate compression ratios >200, duplicate `.dist-info`, missing metadata, `RECORD.jws`/`RECORD.p7s`, duplicate/missing/extra RECORD rows, non-UTF-8 CSV, padded or non-URL-safe base64, non-SHA256 hashes, wrong size/hash, and non-empty hash/size on the RECORD row.

  Add contract cases for missing packages/data, wrong entry points, raw scripts, fixtures, installed legacy paths, caches, multiple native extensions, or a native extension outside `chiq/_bse_solver<suffix>`.

- [ ] **Step 2: Run the tests and confirm red.**

  ```bash
  python3 -m pytest tests/python/non-mpi/packaging/test_wheel_inspector.py -q
  ```

  Expected failure: wheel-validation functions are absent.

- [ ] **Step 3: Implement structural and cryptographic validation.**

  Add:

  ```python
  def validate_wheel(path):
      """Validate ZIP safety, wheel metadata, RECORD, and ChiQ content."""
  ```

  Check ZIP central-directory offsets/ranges and reject overlap before reading entries. Read every member fully to force CRC validation while enforcing logical limits. Parse RECORD as PEP 376/427 UTF-8 CSV. For every non-RECORD member require `sha256=<urlsafe-base64-without-padding>` and exact byte length; decode by restoring padding only internally, then compare `hashlib.sha256(data).digest()`.

  Require exactly one normalized `chiq-<version>.dist-info`; compare wheel filename, METADATA name/version, and dist-info name using `packaging.utils.canonicalize_name` and `parse_wheel_filename`. Require Python `>=3.8`, dependency/extra declarations, and the exact 22 entry points.

  Require `chiq`, `chiq/cli`, `chiq/point_group_data`, `bse`, and `bse_solver`. Allow exactly one native member whose basename begins `_bse_solver` and whose parent is `chiq`; reject every other native suffix. This is the artifact-level allowlist, not a substitute for runtime import.

- [ ] **Step 4: Run tests and commit.**

  ```bash
  python3 -m pytest tests/python/non-mpi/packaging/test_wheel_inspector.py -q
  git add tests/python/non-mpi/packaging/artifact_inspector.py \
    tests/python/non-mpi/packaging/test_wheel_inspector.py
  git diff --cached --check
  git commit -m "test(packaging): enforce wheel archive contracts"
  ```

### Task 5: Orchestrate isolated wheel, sdist, and editable verification

**Files:**

- Create: `tests/python/non-mpi/packaging/runtime_smoke.py`
- Create: `tests/python/non-mpi/packaging/test_runtime_smoke.py`
- Modify: `tests/python/non-mpi/packaging/verify_pip_install.sh`
- Modify: `pyproject.toml`

- [ ] **Step 1: Write runtime and harness-shape tests.**

  Unit-test origin assertions for wheel/sdist venvs and editable mode. Require the exact ABI call:

  ```python
  solver = BSESolver(
      12.0,
      np.array([2], dtype=int),
      np.array([1, 1], dtype=int),
      np.array([1], dtype=int),
  )
  solver.set({
      (0, 0): np.ones((1, 1), dtype=np.complex128),
      (1, 1): np.ones((1, 1), dtype=np.complex128),
  }, "X0_loc")
  ```

  Extend `test_packaging_ci.py` or add shell-source assertions proving three distinct venvs, hash rechecks before install, build-from-extracted-sdist, `--no-cache-dir`, external working directories, and no `cp -R` of the checkout.

- [ ] **Step 2: Run tests and confirm red.**

  ```bash
  python3 -m pytest tests/python/non-mpi/packaging/test_runtime_smoke.py \
    tests/python/non-mpi/packaging/test_packaging_ci.py -q
  ```

  Expected failures: the runtime module is absent and the current shell script reuses one dirty checkout and one venv.

- [ ] **Step 3: Implement runtime checks.**

  `runtime_smoke.py` accepts `--mode wheel|sdist|editable`, `--expected-root`, and `--snapshot-root`. From an unrelated current directory it must:

  - import `chiq`, `chiq._bse_solver`, `bse`, and `bse_solver`;
  - assert wheel/sdist origins are inside the active venv;
  - assert editable pure-Python origins map to the snapshot, the extension does not live under `python/package/chiq`, and the snapshot has no native suffixes;
  - execute the ABI smoke, all shim checks, all canonical help/version calls, and all alias forwarding/warning checks;
  - prove core imports do not load matplotlib or DCore.

- [ ] **Step 4: Rewrite the authoritative pip verifier.**

  Start with `env -i` semantics implemented through an explicit allowlist. Scrub at least `PYTHONPATH`, `PYTHONHOME`, `VIRTUAL_ENV`, `CONDA_PREFIX`, `CMAKE_PREFIX_PATH`, `CMAKE_BUILD_PARALLEL_LEVEL`, `SKBUILD_*`, `PIP_*` except the intentional no-cache setting, and checkout/build path entries. Restore only `HOME`, locale, system tool paths, `TMPDIR`, selected compiler variables (`CC`, `CXX`, `SDKROOT`), `CMAKE_ARGS`, and `EIGEN3_INCLUDE_DIR`; print the resulting named environment without secrets.

  The sequence is mandatory:

  1. Create the Git-backed snapshot.
  2. Create a frontend venv; on Python 3.8 install `pip<25.1` and `build<1.3`, otherwise current supported frontends.
  3. Build an sdist into a retained artifact directory with verbose logs, hash it, validate it, and safely extract it.
  4. Build the wheel only from that extracted gitless/submoduleless sdist, retain/hash/validate it, and capture verbose backend logs.
  5. Create fresh venv A, re-hash and install that exact wheel with `--no-cache-dir`, then run runtime smoke externally.
  6. Create fresh venv B, re-hash and install that exact sdist, then run runtime smoke externally.
  7. Create fresh venv C, install the snapshot editably, then run runtime smoke externally.

  Write `python -m pip list --format=json`, frontend versions, parsed build-system distribution versions, hashes, and snapshot manifest to the retained diagnostics directory. Each mode gets separate build/cache/temp directories and cannot see artifacts from another mode.

- [ ] **Step 5: Add build exclusions without hiding the authoritative extension.**

  Keep the existing `wheel.packages` exactly. Add `sdist.exclude` for package-root native files, caches, and build trees. Add `wheel.exclude` for cache/bytecode content only: do **not** exclude `chiq/_bse_solver*`, because that would also remove the sole CMake-produced extension. The clean snapshot plus wheel allowlist rejects source-tree native contamination. Add a metadata test documenting this distinction.

- [ ] **Step 6: Run the unit suite, then one local endpoint flow.**

  ```bash
  python3 -m pytest tests/python/non-mpi/packaging/test_runtime_smoke.py \
    tests/python/non-mpi/packaging/test_packaging_ci.py \
    tests/python/non-mpi/packaging/test_pyproject_metadata.py -q
  PYTHON="$(command -v python3)" \
    tests/python/non-mpi/packaging/verify_pip_install.sh
  ```

  Expected final line: `ALL PIP INSTALL PATHS VERIFIED`. Inspect retained logs to confirm the wheel build source is the extracted sdist and backend versions were recorded.

- [ ] **Step 7: Commit.**

  ```bash
  git add pyproject.toml tests/python/non-mpi/packaging/runtime_smoke.py \
    tests/python/non-mpi/packaging/test_runtime_smoke.py \
    tests/python/non-mpi/packaging/test_packaging_ci.py \
    tests/python/non-mpi/packaging/test_pyproject_metadata.py \
    tests/python/non-mpi/packaging/verify_pip_install.sh
  git diff --cached --check
  git commit -m "test(packaging): isolate pip artifact verification"
  ```

### Task 6: Make legacy native cleanup narrow and testable

**Files:**

- Create: `cmake/RemoveObsoleteLegacySolver.cmake.in`
- Modify: `CMakeLists.txt`
- Modify: `src/CMakeLists.txt`
- Modify: `python/CMakeLists.txt`
- Create: `tests/python/non-mpi/packaging/test_legacy_cleanup.py`
- Modify: `tests/python/non-mpi/packaging/verify_legacy_install.sh`

- [ ] **Step 1: Write cleanup tests before CMake changes.**

  Configure tiny temporary installs covering absolute and relative `CMAKE_INSTALL_LIBDIR`, configure-time `CMAKE_INSTALL_PREFIX`, `DESTDIR`, a missing target directory, obsolete regular files and symlinks, matching directories/other types, symlink ancestors, and adjacent sentinel files. Require missing target directories to be a successful no-op. Require failures to leave sentinels unchanged.

- [ ] **Step 2: Run tests and confirm red.**

  ```bash
  python3 -m pytest tests/python/non-mpi/packaging/test_legacy_cleanup.py -q
  ```

  Expected failure: the cleanup template and install hook do not exist.

- [ ] **Step 3: Implement install-time cleanup immediately before target install.**

  Configure a script containing the configured full libdir. At install time combine `$ENV{DESTDIR}` with that path, normalize with `get_filename_component(... ABSOLUTE)`, append `bse-python`, and verify every ancestor with `IS_SYMLINK` before globbing. If the exact directory is absent, return successfully. Remove only regular files or symlinks matching:

  ```cmake
  bse_solver*.so  bse_solver*.pyd  bse_solver*.dylib
  chiq/_bse_solver*.so  chiq/_bse_solver*.pyd  chiq/_bse_solver*.dylib
  ```

  Fail on a matching directory or other type and verify every candidate remains inside the intended directory. Register the install script immediately before `install(TARGETS _bse_solver ...)`. Acceptance uses configure-time `-DCMAKE_INSTALL_PREFIX=...` followed by `make install`; do not depend on `cmake --install --prefix` override semantics.

  In `python/CMakeLists.txt`, add `PATTERN` exclusions for `__pycache__`, bytecode, and native suffixes to every `install(DIRECTORY)`. Remove only the exact obsolete `bse` symlink immediately before installing the real directory.

- [ ] **Step 4: Expand the legacy verifier.**

  Use a temporary build directory outside the repository; never delete `build_legacy`. Configure `CHIQ_WHEEL_BUILD=OFF` against the real checkout and populated pybind11 v2.13.6 submodule. Enable and run C++ tests. Install into a temporary configure-time prefix and, from an unrelated directory, assert:

  - `chiq`, `bse`, `bse_solver`, and the sole extension resolve under the install prefix;
  - `bse` is a real directory;
  - all 11 suffixed and 11 suffix-free executables exist and are executable;
  - every executable passes help/version, aliases warn, and the ABI smoke passes;
  - legacy C++ tests and the complete non-MPI suite pass against installed origins.

- [ ] **Step 5: Run verification and commit.**

  ```bash
  python3 -m pytest tests/python/non-mpi/packaging/test_legacy_cleanup.py -q
  tests/python/non-mpi/packaging/verify_legacy_install.sh
  git add CMakeLists.txt cmake/RemoveObsoleteLegacySolver.cmake.in \
    src/CMakeLists.txt python/CMakeLists.txt \
    tests/python/non-mpi/packaging/test_legacy_cleanup.py \
    tests/python/non-mpi/packaging/verify_legacy_install.sh
  git diff --cached --check
  git commit -m "fix(build): clean obsolete legacy solver binaries safely"
  ```

### Task 7: Pin the DCore boundary and Python 3.8-compatible build frontends

**Files:**

- Modify: `pyproject.toml`
- Create: `tests/python/non-mpi/packaging/dcore_boundary_smoke.py`
- Create: `tests/python/non-mpi/packaging/test_dcore_boundary.py`
- Modify: `tests/python/non-mpi/packaging/test_pyproject_metadata.py`
- Modify: `.github/workflows/main.yml`

- [ ] **Step 1: Write metadata and boundary tests.**

  Require these effective declarations:

  ```toml
  requires = [
    "scikit-build-core>=0.10,<0.11; python_version < '3.9'",
    "scikit-build-core>=0.10,<1; python_version >= '3.9'",
    "pybind11>=2.12,<3",
    "cmake>=3.15",
  ]
  dcore = ["dcore==4.2.0", "mpi4py"]
  mpi = ["mpi4py"]
  ```

  Unit-test missing-extra behavior for canonical and deprecated `dcore_chiq` commands. Assert the workflow contains endpoint jobs for 3.8 and 3.13, OpenMPI installation, exact-extra installation, no-cache settings, and retained failure diagnostics.

- [ ] **Step 2: Run tests and confirm red.**

  ```bash
  python3 -m pytest tests/python/non-mpi/packaging/test_dcore_boundary.py \
    tests/python/non-mpi/packaging/test_pyproject_metadata.py \
    tests/python/non-mpi/packaging/test_packaging_ci.py -q
  ```

  Expected failures: DCore is unpinned, mpi4py is missing from the DCore extra, build requirements lack 3.8 markers, and endpoint boundary jobs are absent.

- [ ] **Step 3: Implement the exact DCore import boundary smoke.**

  `dcore_boundary_smoke.py` imports and checks callable/presence for the exact symbols used by `chiq.cli.dcore_chiq` in `dcore._dispatcher`, `dcore.dmft_core`, `dcore.program_options`, `dcore.tools`, `dcore.impurity_solvers`, and `dcore.sumkdft_workers.launcher`, then imports `mpi4py`. Do not assert private signatures.

  Run canonical and alias help/version, then invoke a nonexistent input path that completes all DCore imports and reaches post-import input validation. Assert a controlled nonzero exit without traceback from an import failure. This replaces, rather than supplements, the base design's numerical DCore end-to-end promise.

- [ ] **Step 4: Partition CI by endpoint and failure domain.**

  On both Python 3.8 and 3.13:

  - run the authoritative pip flow exactly once with caches disabled;
  - run a core-only missing-DCore-extra job;
  - run an Ubuntu/OpenMPI DCore-boundary job with `.[dcore]`;
  - upload manifests, artifacts, and verbose logs with `if: failure()`.

  Keep the Python 3.11 smoke job. Use job-level timeouts so archive/unit, pip-artifact, legacy, and DCore failures are distinguishable. Exact dependency drift or Python 3.8 ecosystem failure is intentionally surfaced, not silently relaxed.

- [ ] **Step 5: Run local static tests and commit.**

  ```bash
  python3 -m pytest tests/python/non-mpi/packaging/test_dcore_boundary.py \
    tests/python/non-mpi/packaging/test_pyproject_metadata.py \
    tests/python/non-mpi/packaging/test_packaging_ci.py -q
  git add pyproject.toml .github/workflows/main.yml \
    tests/python/non-mpi/packaging/dcore_boundary_smoke.py \
    tests/python/non-mpi/packaging/test_dcore_boundary.py \
    tests/python/non-mpi/packaging/test_pyproject_metadata.py \
    tests/python/non-mpi/packaging/test_packaging_ci.py
  git diff --cached --check
  git commit -m "ci(packaging): verify endpoint dependency boundaries"
  ```

### Task 8: Document the supported installation and verification matrix

**Files:**

- Modify: `README.md`
- Modify: `doc/install.rst`
- Modify: `tests/python/non-mpi/packaging/test_pyproject_metadata.py`

- [ ] **Step 1: Add a failing documentation assertion.**

  Require user docs to state Python >=3.8, continuous Python 3.8/3.13 verification, wheel/sdist/editable/legacy coverage, the isolated `chiq[dcore]` recommendation with DCore 4.2.0, and `chiq[plot]` for plotting commands.

- [ ] **Step 2: Run and confirm red.**

  ```bash
  python3 -m pytest tests/python/non-mpi/packaging/test_pyproject_metadata.py -q
  ```

  Expected failure: one or more matrix/extra statements are missing.

- [ ] **Step 3: Update docs without promising unsupported targets.**

  State that portable wheel-member validation recognizes `.pyd`, but Windows snapshot no-follow behavior and legacy cleanup are not acceptance targets. Document that unsupported no-follow platforms fail closed. Keep Python 3.8 support explicit despite upstream EOL, and explain that endpoint CI detects dependency drift.

- [ ] **Step 4: Build docs, run tests, and commit.**

  ```bash
  python3 -m pytest tests/python/non-mpi/packaging/test_pyproject_metadata.py -q
  python3 -m sphinx -b html doc /tmp/chiq-packaging-docs
  git add README.md doc tests/python/non-mpi/packaging/test_pyproject_metadata.py
  git diff --cached --check
  git commit -m "docs: describe packaging verification guarantees"
  ```

  Expected docs result: build succeeds; the pre-existing missing `support/faqs` warning may remain, but no new warnings are introduced.

### Task 9: Run the complete acceptance matrix and request final review

**Files:**

- Modify only if a verification failure proves a defect in an earlier task; commit each fix separately with its regression test.

- [ ] **Step 1: Run fast static and adversarial suites.**

  ```bash
  python3 -m pytest tests/python/non-mpi/packaging -q
  python3 -m pyflakes python/package/chiq/cli tests/python/non-mpi/packaging
  git diff --check
  ```

- [ ] **Step 2: Run the complete non-MPI suite.**

  ```bash
  cd build
  PYTHONPATH="$OLDPWD/python/package:$OLDPWD/build/src" \
    PATH="$OLDPWD/python/scripts:$PATH" \
    python3 -m pytest ../tests/python/non-mpi -q
  cd ..
  ```

- [ ] **Step 3: Run local pip and legacy acceptance flows.**

  Export the exact local toolchain flags from `AGENTS.md`, then run:

  ```bash
  tests/python/non-mpi/packaging/verify_pip_install.sh
  tests/python/non-mpi/packaging/verify_legacy_install.sh
  ```

  Expected final lines: `ALL PIP INSTALL PATHS VERIFIED` and `LEGACY INSTALL VERIFIED`. Confirm the exact retained wheel and sdist hashes did not change between validation and installation.

- [ ] **Step 4: Confirm repository hygiene.**

  ```bash
  git status --short
  git log --oneline --decorate -10
  ```

  Expected: only the pre-existing untracked scientific outputs remain; no wheel, sdist, cache, native extension, temporary build directory, or local toolchain setting is tracked.

- [ ] **Step 5: Request two independent code reviews.**

  Use `superpowers:requesting-code-review` and the available AI review cycle. Ask one reviewer to check security/architecture (descriptor walking, archive graph validation, RECORD, cleanup containment) and the other to check workflow/user behavior (CLI, shims, extras, CI, docs). Resolve every must-fix with a failing regression test, rerun the affected acceptance flow, and commit the correction.

- [ ] **Step 6: Run final verification after the last review fix.**

  ```bash
  python3 -m pytest tests/python/non-mpi/packaging -q
  git diff --check
  git status --short
  ```

  Do not claim completion from earlier output; the final report must quote the fresh test counts and both authoritative verifier results.

## Specification coverage audit

| Design requirement | Implemented/tested in |
|---|---|
| Descriptor-safe Git snapshot, races, modes, manifest | Task 2 |
| Gitless/submoduleless independent artifact flow | Task 5 |
| Hostile tar validation and manual extraction | Task 3 |
| Hostile ZIP, metadata, RECORD, exact native allowlist | Task 4 |
| Exact sdist contents including HDF5 | Task 3 |
| Existing three-package mapping and build exclusions | Task 5 |
| Narrow legacy cleanup, DESTDIR, sentinels, no-op | Task 6 |
| Shim object equivalence and independent warnings | Task 1 |
| 11 canonical + 11 alias CLI behavior | Tasks 1 and 5 |
| Plot missing-extra UX | Task 1 |
| DCore 4.2.0/OpenMPI boundary at both endpoints | Task 7 |
| NumPy/pybind11 `BSESolver.set` smoke | Tasks 5 and 6 |
| Wheel/sdist/editable/legacy origin assertions | Tasks 5 and 6 |
| Python 3.8 build markers and frontend constraints | Tasks 5 and 7 |
| CI partitioning, diagnostics, docs | Tasks 7 and 8 |
| Full verification and independent review | Task 9 |
