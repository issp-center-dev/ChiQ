import importlib
import sys
import warnings


def test_chiq_has_underscore_bse_solver():
    m = importlib.import_module("chiq._bse_solver")
    assert hasattr(m, "BSESolver")


def test_top_level_bse_solver_shim_warns_and_forwards():
    sys.modules.pop("bse_solver", None)
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        import bse_solver
    from chiq import _bse_solver

    assert bse_solver.BSESolver is _bse_solver.BSESolver
    assert any(issubclass(x.category, FutureWarning) for x in w)
