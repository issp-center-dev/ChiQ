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


class BoundaryError(RuntimeError):
    """An installed dependency or command violates the tested boundary."""


def _require(condition, message):
    if not condition:
        raise BoundaryError(message)


def _inside(path, root):
    path = Path(path).resolve()
    root = Path(root).resolve()
    try:
        path.relative_to(root)
    except ValueError:
        return False
    return True


def _assert_installed_origin(name, origin, expected_root):
    _require(origin, "%s has no filesystem origin" % name)
    _require(
        _inside(origin, expected_root),
        "%s origin is outside the install environment: %s" % (name, origin),
    )


def _check_import_boundary(expected_root):
    dcore = importlib.import_module("dcore")
    _assert_installed_origin("dcore", getattr(dcore, "__file__", None), expected_root)
    version = importlib.metadata.version("dcore")
    _require(version == "4.2.0", "expected DCore 4.2.0, found %r" % version)

    for module_name, names in REQUIRED_CALLABLES.items():
        module = importlib.import_module(module_name)
        _assert_installed_origin(
            module_name, getattr(module, "__file__", None), expected_root
        )
        for name in names:
            value = getattr(module, name, None)
            _require(callable(value), "%s.%s must be callable" % (module_name, name))

    for module_name, names in REQUIRED_MEMBERS.items():
        module = importlib.import_module(module_name)
        _assert_installed_origin(
            module_name, getattr(module, "__file__", None), expected_root
        )
        for name in names:
            _require(hasattr(module, name), "%s.%s is missing" % (module_name, name))

    mpi4py_version = importlib.metadata.version("mpi4py")
    _require(mpi4py_version, "mpi4py distribution has no version")

    impurity_solvers = importlib.import_module("dcore.impurity_solvers")
    _require(
        callable(impurity_solvers.compute_basis_rot),
        "dcore.impurity_solvers.compute_basis_rot must be callable",
    )
    _require(
        isinstance(impurity_solvers.solver_classes, dict),
        "dcore.impurity_solvers.solver_classes must be a dict",
    )


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
        _require(path, "missing installed command: %s" % command)
        _assert_installed_origin(command, path, expected_root)
        paths[command] = path
    return paths


def _check_help_and_version(commands, cwd, version):
    expected_version = "ChiQ version %s\n" % version
    for command, executable in commands.items():
        for option in ("--help", "--version"):
            result = _run((executable, option), cwd)
            combined = result.stdout + result.stderr
            _require(
                result.returncode == 0,
                "%s %s failed:\n%s" % (command, option, combined),
            )
            _require(
                "Traceback" not in combined,
                "%s %s emitted a traceback:\n%s" % (command, option, combined),
            )
            if option == "--help":
                _require(
                    "usage:" in result.stdout.lower(),
                    "%s --help omitted usage output" % command,
                )
            else:
                _require(
                    result.stdout == expected_version,
                    "%s --version output was %r, expected %r"
                    % (command, result.stdout, expected_version),
                )
            if command.endswith(".py"):
                _require(
                    "deprecated" in result.stderr.lower(),
                    "%s omitted its deprecation warning" % command,
                )
            else:
                _require(
                    result.stderr == "",
                    "%s emitted unexpected stderr: %s" % (command, result.stderr),
                )


def _check_post_import_validation(commands, cwd, mode):
    missing = Path(cwd) / "chiq-dcore-boundary-input-does-not-exist.ini"
    _require(
        not missing.exists(), "boundary fixture unexpectedly exists: %s" % missing
    )
    for command, executable in commands.items():
        result = _run((executable, os.fspath(missing), "--np", "1"), cwd)
        combined = result.stdout + result.stderr
        _require(result.returncode != 0, "%s unexpectedly succeeded" % command)
        _require(
            "Traceback" not in combined,
            "%s input validation emitted a traceback:\n%s" % (command, combined),
        )
        if mode == "missing-extra":
            _require(
                "pip install chiq[dcore]" in result.stderr,
                "%s omitted the missing-extra instruction:\n%s" % (command, combined),
            )
        else:
            _require(
                "does not exist" in result.stderr,
                "%s did not reach input validation:\n%s" % (command, combined),
            )
            _require(
                "pip install chiq[dcore]" not in combined,
                "%s reported the installed DCore extra as missing:\n%s"
                % (command, combined),
            )
            _require(
                "ImportError" not in combined,
                "%s failed at the DCore import boundary:\n%s" % (command, combined),
            )
        if command.endswith(".py"):
            _require(
                "deprecated" in result.stderr.lower(),
                "%s omitted its deprecation warning" % command,
            )


def _main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=("installed", "missing-extra"), required=True)
    parser.add_argument("--expected-root", required=True)
    arguments = parser.parse_args(argv)

    import chiq

    _assert_installed_origin(
        "chiq", getattr(chiq, "__file__", None), arguments.expected_root
    )
    commands = _command_paths(arguments.expected_root)
    cwd = Path.cwd().resolve()
    _check_help_and_version(commands, cwd, chiq.__version__)
    if arguments.mode == "installed":
        _check_import_boundary(arguments.expected_root)
    _check_post_import_validation(commands, cwd, arguments.mode)
    print("DCore boundary smoke OK (%s)" % arguments.mode)
    return 0


def main(argv=None):
    try:
        return _main(argv)
    except (BoundaryError, ImportError, importlib.metadata.PackageNotFoundError) as error:
        print("DCore boundary verification failed: %s" % error, file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
