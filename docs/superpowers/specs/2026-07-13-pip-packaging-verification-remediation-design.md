# Pip-packaging verification remediation design

**Date:** 2026-07-13  
**Status:** Approved; AI review-cycle Phase 0 converged  
**Extends:** `2026-07-12-pip-packaging-design.md`, especially §7

## Goal and scope

Close the remaining packaging-verification and artifact-isolation gaps while preserving:

- Python 3.8 as the minimum supported version;
- the legacy direct-CMake workflow;
- the wheel, sdist, and editable-install workflows;
- all scientific calculations and default CLI behavior.

The only additive CLI behavior is uniform `--version` support. This addendum supersedes any
ambiguous acceptance wording in §7 of the base design.

## 1. Snapshot boundary

The build backend never invokes Git. The authoritative verification harness explicitly requires a
Git checkout and invokes `git ls-files --stage -z` and `git diff HEAD --stat`; Git is not a build
dependency.

The snapshot builder:

- rejects unmerged and non-stage-0 entries;
- skips mode `160000` gitlinks;
- rejects mode `120000` symlinks and unknown modes;
- accepts only `100644` and `100755`, with worktree type and executable mode matching the index;
- uses `lstat` on every ancestor and leaf, rejecting symlink/non-directory ancestors and
  non-regular leaves;
- opens leaves with `O_NOFOLLOW`, then compares `fstat` device, inode, and type with `lstat`,
  failing on races or type changes;
- fails closed when an equivalent no-follow operation is unavailable.

It copies current tracked-worktree bytes, intentionally including unstaged modifications to
tracked files. New files must be staged. It warns about excluded untracked paths and records HEAD,
`git diff HEAD --stat`, and a sorted SHA-256/index-mode/path manifest. The resulting snapshot
contains no `.git`, submodule contents, native binaries, or caches.

## 2. Artifact flow

Python 3.8 and 3.13 each execute the following flow independently. `PYTHONPATH`, `PYTHONHOME`,
checkout/build environment variables, and checkout/build `PATH` entries are removed first.

1. Build, retain, inspect, and hash an sdist from the snapshot.
2. Safely extract it into a new gitless, submodule-less directory.
3. Build, retain, inspect, and hash a wheel from that sdist. Wheel mode uses the existing
   `find_package(pybind11)` path and build-system dependency.
4. Recheck the hash and install that exact wheel with `--no-cache-dir` into fresh venv A.
5. Recheck and install that exact sdist into fresh venv B.
6. Install the snapshot editably into fresh venv C.
7. Run every acceptance check from outside all source and build trees.

No environment or native artifact is reused between modes or interpreters.

## 3. Archive validation and extraction

Before parsing, reject an archive larger than 512 MiB. Logical limits are 50,000 members,
256 MiB per file, and 2 GiB total uncompressed data.

Archive paths are POSIX-only. Reject NULs, backslashes, absolute/drive/UNC paths, repeated
separators, empty/`.`/`..` components, non-canonical paths, trailing aliases, duplicate identities,
file/directory conflicts, and ancestor-file conflicts.

### Tar sdists

- Accept only regular files and directories; reject links, devices, FIFOs, sparse files, and
  special members.
- Require one shared normalized name/version top-level component; an explicit root entry is not
  required.
- Enforce aggregate compression ratio ≤200; tar does not use a per-member ratio.
- Fully resolve both extraction root and destination before containment checks.
- Create a new root with mode `0700`, validate the complete graph before writing, and never call
  `extractall`.
- Create files exclusively. Apply `0755` only to source owner-executable files and `0644`
  otherwise; directories become `0755`. Never restore ownership or timestamps.

### Wheel ZIP files

Wheel ZIPs are validated but not manually extracted. Reject encrypted, multi-disk, overlapping,
out-of-range, or CRC-invalid members. Enforce per-entry and aggregate compression ratios ≤200.
Forbid `RECORD.jws` and `RECORD.p7s`.

Unit tests cover traversal, links, special members, archive bombs, aliases, conflicts, permissions,
and macOS `/tmp` resolution.

## 4. Wheel and sdist contracts

The wheel must have exactly one `.dist-info` directory with `METADATA`, `entry_points.txt`, and
`RECORD`. Filename and metadata name/version are compared with standard packaging normalization.
There must be exactly one `RECORD` row per file. Only `RECORD` itself may omit hash and size; every
other row requires a valid SHA-256 and exact size.

The wheel must contain:

- `chiq`, including `cli`, `point_group_data`, and exactly one `chiq._bse_solver` extension;
- `bse` and `bse_solver`;
- metadata declaring Python ≥3.8, all dependencies/extras, 11 canonical entry points, and their
  11 deprecated `.py` aliases.

It must not contain tests/HDF5 fixtures, raw `python/scripts`, `lib/bse-python`, `chiqvars.sh`,
caches, build/source trees, or duplicate/non-`chiq` solver extensions. Native suffix matching never
replaces an actual import and ABI smoke test.

The sdist must contain `pyproject.toml`, `LICENSE`, `README.md`, all required CMake files, C++
sources/headers, CMake helpers, the three Python package trees, CLI and point-group modules, and
test HDF5 fixtures. It must exclude native binaries, caches, build trees, and installed layouts.

## 5. Package mapping and legacy cleanup

Keep the existing scikit-build-core mapping:

```toml
wheel.packages = [
  "python/package/chiq",
  "python/package/bse",
  "python/package/bse_solver",
]
```

Add scikit-build and `install(DIRECTORY)` exclusions for caches and native artifacts.
`src/CMakeLists.txt` remains the sole native-install authority.

Immediately before legacy `install(TARGETS)`, an install-time CMake script:

- combines install-time `DESTDIR` with `CMAKE_INSTALL_FULL_LIBDIR`, normalizes the result with
  `get_filename_component(... ABSOLUTE)`, and appends `bse-python`;
- removes only `bse_solver*.{so,pyd,dylib}` in that exact directory;
- removes only `_bse_solver*.{so,pyd,dylib}` in its exact `chiq/` child;
- removes only matching regular files or symlinks, fails on matching directories/other types,
  and verifies containment;
- installs the current target afterward.

The exact old `bse` symlink is removed immediately before installing the real `bse` directory.
Tests cover relative/absolute install libdirs, `DESTDIR`, obsolete regular/symlink files, and
neighboring sentinels that must survive.

## 6. Compatibility shims

The compatibility manifest covers:

- `bse_toml`, `g2scl_core`, `h5bse`, `index_pair`, `matrix_dict`, `mpi`, `point_group`,
  `sumk_dft_chi`, and `tools`;
- point-group modules `C1`, `C2`, `D3`, `D4`, `D6`, `O`, `Oh`, and `base`;
- `bse_solver.BSESolver`.

Forwarding modules are distinct module objects. Tests assert only named public exported-object
equivalence and behavior. In a fresh process, the first `bse` or nested `bse` import emits one
`bse` warning total. `bse_solver` independently emits its own one warning, including after `bse`.

## 7. CLI and optional extras

All 11 canonical argparse commands gain a shared `--version` option sourced from
`chiq.__version__`. It prints exactly `ChiQ version <version>` and exits successfully. Target-module
and transitive imports, `--help`, and `--version` must work in a core-only environment. Help uses
the invoked `argv[0]` rather than a hard-coded deprecated name.

Existing non-MPI tests plus a valid existing invocation for every changed parser protect argument
semantics. Plot commands retain argparse-first invalid-syntax behavior; syntactically valid minimal
fixtures reach the matplotlib gate and, without the extra, produce actionable stderr naming
`pip install chiq[plot]` with no traceback.

Wheel `.py` aliases are executable `[project.scripts]` entry points. Legacy CMake separately
installs executable wrapper files. Every alias warns on stderr and forwards unchanged arguments.

## 8. DCore boundary

Because ChiQ imports private DCore APIs, the published extra pins `dcore==4.2.0` and includes
`mpi4py`; the independent `mpi` extra remains `mpi4py`. Documentation recommends an isolated
environment. DCore upgrades require a reviewed compatibility change.

Ubuntu/OpenMPI CI on Python 3.8 and 3.13 verifies the exact imported symbols in:

- `dcore._dispatcher`, `dcore.dmft_core`, `dcore.program_options`, and `dcore.tools`;
- `dcore.impurity_solvers` and `dcore.sumkdft_workers.launcher`;
- `mpi4py`.

The promise is import/callable presence, not private-function signatures or scientific execution.
The job runs canonical and deprecated `dcore_chiq` help/version plus a nonexistent-input invocation
that reaches post-import validation. Core-only environments test the corresponding missing-extra
message on both endpoint Pythons.

## 9. ABI, editable, and legacy checks

The native ABI smoke constructs:

```python
BSESolver(
    12.0,
    np.array([2], dtype=int),
    np.array([1, 1], dtype=int),
    np.array([1], dtype=int),
)
```

It then sets `X0_loc` with keys `(0, 0)` and `(1, 1)`, each holding a `complex128` 1×1 array of
ones. This crosses the NumPy/pybind11 boundary that previously segfaulted.

Wheel/sdist module origins must lie in their fresh venv. In editable mode, pure Python maps to the
snapshot, the native extension must not lie under `python/package/chiq`, the snapshot contains no
native binary, and the ABI smoke must pass.

Legacy verification runs only from the real checkout with the populated pybind11 v2.13.6
submodule and `CHIQ_WHEEL_BUILD=OFF`. From an unrelated directory it verifies installed origins,
the real `bse` directory, one native extension, all 11 suffixed and 11 suffix-free executables,
help/version, ABI smoke, C++ tests, and non-MPI tests.

## 10. Python 3.8, CI, and documentation

Build-system requirements use environment markers:

- Python <3.9: `scikit-build-core>=0.10,<0.11`;
- Python ≥3.9: `scikit-build-core>=0.10,<1`;
- all versions: `pybind11>=2.12,<3` and `cmake>=3.15`.

The Python 3.8 harness uses `pip<25.1` and `build<1.3`. Retained manifests record exact frontend
and backend versions. The 3.8/3.13 matrix executes the authoritative flow once per endpoint with
no cache and uploads artifacts/logs on failure. The 3.11 smoke job remains. Documentation states
Python ≥3.8 support and continuous 3.8/3.13 validation.

Windows snapshot and legacy cleanup are not acceptance targets; unsupported no-follow behavior
fails closed. Portable wheel-member checks still recognize `.pyd`, and Windows locked-file cleanup
never broadens deletion.

## 11. TDD and scope controls

Write failing tests first for snapshot no-follow/races/modes, archive security and limits,
permissions, `RECORD`, metadata and 22 scripts, exact cleanup/order/containment/sentinels, shim
warnings, CLI behavior, optional-dependency UX, module origins, and CI endpoints. Then implement the
minimum changes and run full verification.

Do not change numerical behavior, defaults, TOML handling, or add local macOS flags to
`pyproject.toml`. Do not commit generated test outputs.

## AI review-cycle record

Phase 0 ran Codex (technical/architecture) and Antigravity (workflow/user lens). After five design
revisions:

- both reviewers reported zero `must_fix_design` findings;
- all remaining questions and risks received explicit dispositions;
- `evaluate_design_convergence` returned `converged: true` with no unresolved blockers.

