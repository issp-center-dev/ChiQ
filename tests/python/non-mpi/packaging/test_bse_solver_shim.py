import importlib
import sys
import warnings

import pytest


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
