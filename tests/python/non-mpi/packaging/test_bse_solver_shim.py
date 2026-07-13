import importlib
import json
import os
import subprocess
import sys
import warnings

import pytest


ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(__file__)))))


def test_chiq_has_underscore_bse_solver():
    try:
        m = importlib.import_module("chiq._bse_solver")
    except ModuleNotFoundError as exc:
        if exc.name != "chiq._bse_solver":
            raise
        pytest.skip("uninstalled legacy build exposes top-level _bse_solver")
    assert hasattr(m, "BSESolver")


def test_top_level_bse_solver_shim_warns_and_forwards():
    sys.modules.pop("bse_solver", None)
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        import bse_solver
    try:
        expected = importlib.import_module("chiq._bse_solver")
    except ModuleNotFoundError as exc:
        if exc.name != "chiq._bse_solver":
            raise
        expected = importlib.import_module("_bse_solver")

    assert bse_solver.BSESolver is expected.BSESolver
    assert any(issubclass(x.category, FutureWarning) for x in w)


def test_bse_solver_warning_is_independent_from_bse_in_fresh_process():
    env = dict(os.environ)
    env["PYTHONPATH"] = os.pathsep.join((
        os.path.join(ROOT, "python", "package"),
        os.path.join(ROOT, "build", "src"),
    ))
    code = """
import importlib
import json
import warnings
with warnings.catch_warnings(record=True) as caught:
    warnings.simplefilter('always')
    importlib.import_module('bse')
    bse_solver = importlib.import_module('bse_solver')
    importlib.reload(bse_solver)
print(json.dumps([(w.category.__name__, str(w.message)) for w in caught]))
"""
    result = subprocess.run(
        [sys.executable, "-c", code],
        cwd=os.path.dirname(ROOT),
        env=env,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    caught = json.loads(result.stdout)
    assert [category for category, message in caught] == ["FutureWarning", "FutureWarning"]
    assert "'bse' package is deprecated" in caught[0][1]
    assert "top-level 'bse_solver' module is deprecated" in caught[1][1]
