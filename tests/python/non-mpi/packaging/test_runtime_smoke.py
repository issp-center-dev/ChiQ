import importlib.util
import os
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace

import numpy as np
import pytest


MODULE_PATH = Path(__file__).with_name("runtime_smoke.py")


def _load_runtime_smoke():
    spec = importlib.util.spec_from_file_location("runtime_smoke", str(MODULE_PATH))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_runtime_smoke_module_exists():
    assert MODULE_PATH.is_file()


def test_runtime_smoke_imports_in_fresh_interpreter():
    env = dict(os.environ)
    env.pop("PYTHONPATH", None)
    code = """
import importlib
import numpy
import sys

sys.modules.pop("importlib.util", None)
if hasattr(importlib, "util"):
    delattr(importlib, "util")
path = sys.argv[1]
namespace = {"__file__": path, "__name__": "runtime_smoke_import_test"}
with open(path, "rb") as stream:
    exec(compile(stream.read(), path, "exec"), namespace)
print("runtime smoke import OK")
"""
    result = subprocess.run(
        [sys.executable, "-c", code, str(MODULE_PATH)],
        cwd=str(MODULE_PATH.parent),
        env=env,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, result.stderr
    assert result.stdout == "runtime smoke import OK\n"
    assert "Traceback" not in result.stdout + result.stderr


def test_wheel_origin_rejects_checkout_module(tmp_path):
    smoke = _load_runtime_smoke()
    venv = tmp_path / "venv"
    checkout = tmp_path / "checkout"
    module = SimpleNamespace(__file__=str(checkout / "python/package/chiq/__init__.py"))

    with pytest.raises(AssertionError, match="outside expected root"):
        smoke.assert_module_origins("wheel", venv, checkout, {"chiq": module})


def test_wheel_origin_accepts_only_venv_modules(tmp_path):
    smoke = _load_runtime_smoke()
    venv = tmp_path / "venv"
    module = SimpleNamespace(__file__=str(venv / "lib/python/site-packages/chiq/__init__.py"))
    smoke.assert_module_origins("sdist", venv, tmp_path / "snapshot", {"chiq": module})


def test_editable_origins_split_python_and_native(tmp_path):
    smoke = _load_runtime_smoke()
    snapshot = tmp_path / "snapshot"
    venv = tmp_path / "venv"
    modules = {
        "chiq": SimpleNamespace(
            __file__=str(snapshot / "python/package/chiq/__init__.py")
        ),
        "bse": SimpleNamespace(__file__=str(snapshot / "python/package/bse/__init__.py")),
        "bse_solver": SimpleNamespace(
            __file__=str(snapshot / "python/package/bse_solver/__init__.py")
        ),
        "chiq._bse_solver": SimpleNamespace(
            __file__=str(venv / "lib/python/site-packages/chiq/_bse_solver.so")
        ),
    }
    smoke.assert_module_origins("editable", venv, snapshot, modules)


def test_editable_rejects_native_module_in_snapshot(tmp_path):
    smoke = _load_runtime_smoke()
    snapshot = tmp_path / "snapshot"
    modules = {
        "chiq": SimpleNamespace(
            __file__=str(snapshot / "python/package/chiq/__init__.py")
        ),
        "chiq._bse_solver": SimpleNamespace(
            __file__=str(snapshot / "python/package/chiq/_bse_solver.so")
        ),
    }
    with pytest.raises(AssertionError, match="native extension is inside editable source"):
        smoke.assert_module_origins("editable", tmp_path / "venv", snapshot, modules)


def test_abi_smoke_uses_exact_constructor_and_set_shape():
    smoke = _load_runtime_smoke()
    calls = []

    class FakeSolver:
        def __init__(self, *args):
            calls.append(("init", args))

        def set(self, *args):
            calls.append(("set", args))

    smoke.run_abi_smoke(FakeSolver)

    _, constructor = calls[0]
    assert constructor[0] == 12.0
    np.testing.assert_array_equal(constructor[1], np.array([2], dtype=int))
    np.testing.assert_array_equal(constructor[2], np.array([1, 1], dtype=int))
    np.testing.assert_array_equal(constructor[3], np.array([1], dtype=int))
    assert all(array.dtype == np.dtype(int) for array in constructor[1:])
    _, set_args = calls[1]
    assert set(set_args[0]) == {(0, 0), (1, 1)}
    assert all(value.shape == (1, 1) for value in set_args[0].values())
    assert all(value.dtype == np.complex128 for value in set_args[0].values())
    assert set_args[1] == "X0_loc"


def test_cli_check_requires_exact_versions_help_and_alias_warning():
    smoke = _load_runtime_smoke()
    seen = []

    def run(command, cwd):
        seen.append((tuple(command), Path(cwd)))
        name, argument = command
        stderr = ""
        if name.endswith(".py"):
            stderr = "warning: deprecated; use '%s' instead.\n" % name[:-3]
        if argument == "--help":
            return SimpleNamespace(
                returncode=0, stdout="usage: %s\n" % name, stderr=stderr
            )
        return SimpleNamespace(returncode=0, stdout="ChiQ version 1.2.3\n", stderr=stderr)

    smoke.check_console_commands(
        ("one", "two"), "1.2.3", Path("/external"), run=run
    )

    assert seen == [
        (("one", "--help"), Path("/external")),
        (("one.py", "--help"), Path("/external")),
        (("one", "--version"), Path("/external")),
        (("one.py", "--version"), Path("/external")),
        (("two", "--help"), Path("/external")),
        (("two.py", "--help"), Path("/external")),
        (("two", "--version"), Path("/external")),
        (("two.py", "--version"), Path("/external")),
    ]


def test_cli_check_rejects_tracebacks_and_wrong_version():
    smoke = _load_runtime_smoke()

    def run(command, cwd):
        del cwd
        if command[1] == "--version":
            return SimpleNamespace(returncode=0, stdout="wrong\n", stderr="Traceback")
        return SimpleNamespace(returncode=0, stdout="usage: x\n", stderr="")

    with pytest.raises(AssertionError):
        smoke.check_console_commands(("one",), "1.2.3", Path("/external"), run=run)


def test_dcore_dependent_shim_is_located_without_execution(tmp_path):
    smoke = _load_runtime_smoke()
    origins = {
        "chiq.sumk_dft_chi": tmp_path / "chiq/sumk_dft_chi.py",
        "bse.sumk_dft_chi": tmp_path / "bse/sumk_dft_chi.py",
    }
    seen = []

    def find_spec(name):
        seen.append(name)
        return SimpleNamespace(origin=str(origins[name]))

    modules = smoke.locate_dcore_dependent_shims(find_spec=find_spec)

    assert smoke.DCORE_DEPENDENT_SHIMS == ("sumk_dft_chi",)
    assert seen == ["chiq.sumk_dft_chi", "bse.sumk_dft_chi"]
    assert {name: Path(spec.origin) for name, spec in modules.items()} == origins
