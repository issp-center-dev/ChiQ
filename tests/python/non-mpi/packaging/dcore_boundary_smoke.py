#!/usr/bin/env python3
"""Verify ChiQ's installed DCore 4.2.0 import and command boundary."""

import argparse
import importlib
import importlib.metadata
import os
from pathlib import Path
import shutil
import subprocess
import sys


REQUIRED_CALLABLES = {
    "dcore._dispatcher": ("HDFArchive", "dyson"),
    "dcore.dmft_core": ("DMFTCoreSolver",),
    "dcore.program_options": (
        "create_parser", "parse_parameters", "print_parameters",
        "delete_parameters", "OptionStatus",
    ),
    "dcore.tools": (
        "raise_if_mpi_imported", "make_empty_dir",
        "launch_mpi_subprocesses", "float_to_complex_array",
    ),
    "dcore.sumkdft_workers.launcher": ("run_sumkdft",),
}

REQUIRED_MEMBERS = {
    "dcore.impurity_solvers": ("solver_classes", "compute_basis_rot"),
    "mpi4py": ("MPI",),
}


def _inside(path, root):
    path = Path(path).resolve()
    root = Path(root).resolve()
    try:
        path.relative_to(root)
    except ValueError:
        return False
    return True


def _assert_installed_origin(name, origin, expected_root):
    assert origin, "%s has no filesystem origin" % name
    assert _inside(origin, expected_root), (
        "%s origin is outside the install environment: %s" % (name, origin)
    )


def _check_import_boundary(expected_root):
    dcore = importlib.import_module("dcore")
    _assert_installed_origin("dcore", dcore.__file__, expected_root)
    version = importlib.metadata.version("dcore")
    assert version == "4.2.0", "expected DCore 4.2.0, found %r" % version

    for module_name, names in REQUIRED_CALLABLES.items():
        module = importlib.import_module(module_name)
        _assert_installed_origin(module_name, module.__file__, expected_root)
        for name in names:
            value = getattr(module, name, None)
            assert callable(value), "%s.%s must be callable" % (module_name, name)

    for module_name, names in REQUIRED_MEMBERS.items():
        module = importlib.import_module(module_name)
        _assert_installed_origin(module_name, module.__file__, expected_root)
        for name in names:
            assert hasattr(module, name), "%s.%s is missing" % (module_name, name)

    mpi4py_version = importlib.metadata.version("mpi4py")
    assert mpi4py_version, "mpi4py distribution has no version"

    impurity_solvers = importlib.import_module("dcore.impurity_solvers")
    assert callable(impurity_solvers.compute_basis_rot)
    assert isinstance(impurity_solvers.solver_classes, dict)


def _run(command, cwd):
    return subprocess.run(
        command,
        cwd=os.fspath(cwd),
        capture_output=True,
        text=True,
    )


def _command_paths(expected_root):
    paths = {}
    for command in ("dcore_chiq", "dcore_chiq.py"):
        path = shutil.which(command)
        assert path, "missing installed command: %s" % command
        _assert_installed_origin(command, path, expected_root)
        paths[command] = path
    return paths


def _check_help_and_version(commands, cwd, version):
    expected_version = "ChiQ version %s\n" % version
    for command, executable in commands.items():
        for option in ("--help", "--version"):
            result = _run((executable, option), cwd)
            combined = result.stdout + result.stderr
            assert result.returncode == 0, "%s %s failed:\n%s" % (
                command, option, combined
            )
            assert "Traceback" not in combined
            if option == "--help":
                assert "usage:" in result.stdout.lower()
            else:
                assert result.stdout == expected_version
            if command.endswith(".py"):
                assert "deprecated" in result.stderr.lower()
            else:
                assert result.stderr == ""


def _check_post_import_validation(commands, cwd, mode):
    missing = Path(cwd) / "chiq-dcore-boundary-input-does-not-exist.ini"
    assert not missing.exists(), "boundary fixture unexpectedly exists: %s" % missing
    for command, executable in commands.items():
        result = _run((executable, os.fspath(missing), "--np", "1"), cwd)
        combined = result.stdout + result.stderr
        assert result.returncode != 0, "%s unexpectedly succeeded" % command
        assert "Traceback" not in combined, combined
        if mode == "missing-extra":
            assert "pip install chiq[dcore]" in result.stderr, combined
        else:
            assert "does not exist" in result.stderr, combined
            assert "pip install chiq[dcore]" not in combined, combined
            assert "ImportError" not in combined, combined
        if command.endswith(".py"):
            assert "deprecated" in result.stderr.lower()


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=("installed", "missing-extra"), required=True)
    parser.add_argument("--expected-root", required=True)
    arguments = parser.parse_args(argv)

    import chiq

    _assert_installed_origin("chiq", chiq.__file__, arguments.expected_root)
    commands = _command_paths(arguments.expected_root)
    cwd = Path.cwd().resolve()
    _check_help_and_version(commands, cwd, chiq.__version__)
    if arguments.mode == "installed":
        _check_import_boundary(arguments.expected_root)
    _check_post_import_validation(commands, cwd, arguments.mode)
    print("DCore boundary smoke OK (%s)" % arguments.mode)
    return 0


if __name__ == "__main__":
    sys.exit(main())
