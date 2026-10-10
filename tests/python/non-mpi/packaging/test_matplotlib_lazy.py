import builtins
import importlib
import sys
import pytest

def test_g2scl_core_imports_without_matplotlib(monkeypatch):
    # Simulate matplotlib being absent: block its import, drop cached modules.
    for name in list(sys.modules):
        if name == "matplotlib" or name.startswith("matplotlib.") or name == "chiq.g2scl_core":
            monkeypatch.delitem(sys.modules, name, raising=False)
    real_import = builtins.__import__
    def blocked_import(name, *args, **kwargs):
        if name == "matplotlib" or name.startswith("matplotlib"):
            raise ImportError("matplotlib is blocked for this test")
        return real_import(name, *args, **kwargs)
    monkeypatch.setattr(builtins, "__import__", blocked_import)
    # Importing the module must NOT import matplotlib at module scope.
    importlib.import_module("chiq.g2scl_core")
