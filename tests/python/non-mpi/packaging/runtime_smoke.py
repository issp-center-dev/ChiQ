#!/usr/bin/env python3
"""Runtime verification for an installed ChiQ artifact."""

import argparse
import importlib
import importlib.util
import os
from pathlib import Path
import shutil
import subprocess
import sys
import warnings

import numpy as np

from contracts import CANONICAL_COMMANDS, NATIVE_SUFFIXES, POINT_GROUP_MODULES, SHIM_MODULES

DCORE_DEPENDENT_SHIMS = ("sumk_dft_chi",)


def _inside(path, root):
    path = Path(path).resolve()
    root = Path(root).resolve()
    try:
        path.relative_to(root)
    except ValueError:
        return False
    return True


def _module_path(name, module):
    origin = getattr(module, "__file__", None) or getattr(module, "origin", None)
    assert origin, "%s has no filesystem origin" % name
    return Path(origin).resolve()


def locate_dcore_dependent_shims(find_spec=importlib.util.find_spec):
    """Locate installed DCore shims without executing their module bodies.

    Task 7 owns the pinned DCore dependency and executable import boundary.
    Task 5 proves these files are present and correctly placed on core+[mpi].
    """
    located = {}
    for name in DCORE_DEPENDENT_SHIMS:
        for package in ("chiq", "bse"):
            qualified = "%s.%s" % (package, name)
            assert qualified not in sys.modules, "%s was executed unexpectedly" % qualified
            spec = find_spec(qualified)
            assert spec is not None and spec.origin, "missing installed shim: %s" % qualified
            assert qualified not in sys.modules, "%s spec lookup executed module" % qualified
            located[qualified] = spec
    assert not any(name == "dcore" or name.startswith("dcore.") for name in sys.modules)
    return located


def assert_module_origins(mode, expected_root, snapshot_root, modules):
    """Enforce installed and editable origin boundaries."""
    expected_root = Path(expected_root).resolve()
    snapshot_root = Path(snapshot_root).resolve() if snapshot_root else None
    if mode in ("wheel", "sdist"):
        for name, module in modules.items():
            path = _module_path(name, module)
            assert _inside(path, expected_root), (
                "%s origin is outside expected root: %s" % (name, path)
            )
            if snapshot_root is not None:
                assert not _inside(path, snapshot_root), (
                    "%s origin resolves into snapshot: %s" % (name, path)
                )
        return

    assert mode == "editable", "unknown mode: %s" % mode
    assert snapshot_root is not None, "editable mode requires snapshot root"
    source_packages = snapshot_root / "python" / "package"
    native_source = source_packages / "chiq"
    for name, module in modules.items():
        path = _module_path(name, module)
        if name == "chiq._bse_solver":
            assert not _inside(path, native_source), (
                "native extension is inside editable source: %s" % path
            )
            assert _inside(path, expected_root), (
                "editable native extension is outside expected root: %s" % path
            )
        else:
            assert _inside(path, source_packages), (
                "%s editable origin is outside snapshot: %s" % (name, path)
            )


def run_abi_smoke(solver_class):
    solver = solver_class(
        12.0,
        np.array([2], dtype=int),
        np.array([1, 1], dtype=int),
        np.array([1], dtype=int),
    )
    solver.set({
        (0, 0): np.ones((1, 1), dtype=np.complex128),
        (1, 1): np.ones((1, 1), dtype=np.complex128),
    }, "X0_loc")


def _run_command(command, cwd):
    return subprocess.run(
        command,
        cwd=os.fspath(cwd),
        capture_output=True,
        text=True,
    )


def _assert_clean_result(result, command):
    combined = result.stdout + result.stderr
    assert result.returncode == 0, "%s failed:\n%s" % (" ".join(command), combined)
    assert "Traceback" not in combined, "%s emitted a traceback" % " ".join(command)


def check_console_commands(commands, version, external_cwd, run=_run_command):
    expected_version = "ChiQ version %s\n" % version
    for name in commands:
        for argument in ("--help", "--version"):
            canonical = [name, argument]
            result = run(canonical, external_cwd)
            _assert_clean_result(result, canonical)
            if argument == "--help":
                assert "usage:" in result.stdout.lower(), result.stdout
            else:
                assert result.stdout == expected_version, result.stdout
            assert result.stderr == "", result.stderr

            alias = [name + ".py", argument]
            result = run(alias, external_cwd)
            _assert_clean_result(result, alias)
            if argument == "--help":
                assert "usage:" in result.stdout.lower(), result.stdout
            else:
                assert result.stdout == expected_version, result.stdout
            lowered = result.stderr.lower()
            assert "deprecated" in lowered and name in result.stderr, result.stderr


def _import_and_verify_shims():
    modules = {}
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        modules["chiq"] = importlib.import_module("chiq")
        modules["chiq._bse_solver"] = importlib.import_module("chiq._bse_solver")
        modules["bse"] = importlib.import_module("bse")
        modules["bse_solver"] = importlib.import_module("bse_solver")

    future = [item for item in caught if issubclass(item.category, FutureWarning)]
    assert len(future) == 2, "expected exactly two shim FutureWarnings"
    messages = [str(item.message) for item in future]
    assert "'bse' package is deprecated" in messages[0]
    assert "top-level 'bse_solver' module is deprecated" in messages[1]

    native = modules["chiq._bse_solver"]
    solver_shim = modules["bse_solver"]
    assert solver_shim.__all__ == ["BSESolver"]
    assert solver_shim.BSESolver is native.BSESolver

    for name in SHIM_MODULES:
        if name in DCORE_DEPENDENT_SHIMS:
            continue
        canonical = importlib.import_module("chiq.%s" % name)
        shim = importlib.import_module("bse.%s" % name)
        modules["chiq.%s" % name] = canonical
        modules["bse.%s" % name] = shim
        public = [item for item in vars(shim) if not item.startswith("_")]
        assert public, "bse.%s has no public contract" % name
        for item in public:
            assert hasattr(canonical, item), "chiq.%s lacks %s" % (name, item)
            assert getattr(shim, item) is getattr(canonical, item)
    modules.update(locate_dcore_dependent_shims())
    for name in POINT_GROUP_MODULES:
        canonical_name = "chiq.point_group_data.%s" % name
        shim_name = "bse.point_group_data.%s" % name
        canonical = importlib.import_module(canonical_name)
        shim = importlib.import_module(shim_name)
        modules[canonical_name] = canonical
        modules[shim_name] = shim
        public = [item for item in vars(shim) if not item.startswith("_")]
        assert public, "%s has no public contract" % shim_name
        for item in public:
            assert hasattr(canonical, item), "%s lacks %s" % (canonical_name, item)
            assert getattr(shim, item) is getattr(canonical, item)
    return modules


def _assert_no_optional_imports():
    assert "matplotlib" not in sys.modules
    assert not any(name == "dcore" or name.startswith("dcore.") for name in sys.modules)


def _assert_snapshot_has_no_native(snapshot_root):
    package = Path(snapshot_root) / "python" / "package"
    contaminants = [
        path for path in package.rglob("*")
        if path.is_file() and path.name.lower().endswith(NATIVE_SUFFIXES)
    ]
    assert not contaminants, "editable snapshot contains native files: %s" % contaminants


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", required=True, choices=("wheel", "sdist", "editable"))
    parser.add_argument("--expected-root", required=True)
    parser.add_argument("--snapshot-root")
    arguments = parser.parse_args(argv)
    if arguments.mode == "editable" and not arguments.snapshot_root:
        parser.error("editable mode requires --snapshot-root")

    modules = _import_and_verify_shims()
    _assert_no_optional_imports()
    assert_module_origins(
        arguments.mode, arguments.expected_root, arguments.snapshot_root, modules
    )
    if arguments.mode == "editable":
        _assert_snapshot_has_no_native(arguments.snapshot_root)
    run_abi_smoke(modules["chiq._bse_solver"].BSESolver)

    external_cwd = Path.cwd().resolve()
    expected_root = Path(arguments.expected_root).resolve()
    for name in CANONICAL_COMMANDS:
        for executable in (name, name + ".py"):
            location = shutil.which(executable)
            assert location, "missing console command: %s" % executable
            assert _inside(location, expected_root), (
                "console command is outside expected root: %s" % location
            )
    check_console_commands(
        CANONICAL_COMMANDS, modules["chiq"].__version__, external_cwd
    )
    print("%s runtime smoke OK" % arguments.mode)
    return 0


if __name__ == "__main__":
    sys.exit(main())
