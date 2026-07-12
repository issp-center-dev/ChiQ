import builtins
import importlib
import sys
import warnings
import pytest

def _fresh(monkeypatch, *drop_prefixes):
    for name in list(sys.modules):
        if any(name == p or name.startswith(p + ".") or name == p for p in drop_prefixes):
            monkeypatch.delitem(sys.modules, name, raising=False)

def test_import_bse_warns_and_pulls_no_optional_deps(monkeypatch):
    _fresh(monkeypatch, "bse")
    real_import = builtins.__import__
    def blocked(name, *a, **k):
        if name.split(".")[0] in ("mpi4py", "dcore", "matplotlib"):
            raise ImportError(f"{name} blocked")
        return real_import(name, *a, **k)
    monkeypatch.setattr(builtins, "__import__", blocked)
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        importlib.import_module("bse")            # must not import optional deps
    assert any(issubclass(x.category, FutureWarning) for x in w)

def test_from_bse_import_public_submodule():
    from bse import h5bse
    import chiq.h5bse
    assert h5bse.h5BSE is chiq.h5bse.h5BSE       # class identity preserved

def test_nested_point_group_data():
    import bse.point_group_data.C1
    import chiq.point_group_data.C1
    # a public symbol is the same object
    assert dir(bse.point_group_data.C1)

def test_all_allowlisted_modules_importable():
    for mod in ["h5bse", "bse_toml", "matrix_dict", "point_group", "tools", "index_pair", "mpi"]:
        importlib.import_module(f"bse.{mod}")
    for x in ["base", "C1", "C2", "D3", "D4", "D6", "O", "Oh"]:
        importlib.import_module(f"bse.point_group_data.{x}")
