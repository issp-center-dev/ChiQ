import builtins
import importlib
import json
import os
import subprocess
import sys
import warnings

from contracts import POINT_GROUP_MODULES, SHIM_MODULES


ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(__file__)))))


def _fresh_python(code):
    env = dict(os.environ)
    env["PYTHONPATH"] = os.pathsep.join((
        os.path.join(ROOT, "python", "package"),
        os.path.join(ROOT, "build", "src"),
    ))
    result = subprocess.run(
        [sys.executable, "-c", code],
        cwd=os.path.dirname(ROOT),
        env=env,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    assert "Traceback" not in result.stderr
    return result

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
    # a public symbol forwarded via `import *` must be the SAME object
    public = [n for n in dir(chiq.point_group_data.C1) if not n.startswith("_")]
    assert public, "chiq.point_group_data.C1 has no public names to compare"
    name = public[0]
    assert getattr(bse.point_group_data.C1, name) is getattr(chiq.point_group_data.C1, name)


def test_nested_imports_emit_exactly_one_bse_warning_in_fresh_process():
    result = _fresh_python("""
import importlib
import json
import warnings
with warnings.catch_warnings(record=True) as caught:
    warnings.simplefilter('always')
    bse = importlib.import_module('bse')
    importlib.import_module('bse.tools')
    importlib.import_module('bse.point_group_data.C1')
    importlib.reload(bse)
print(json.dumps([w.category.__name__ for w in caught]))
""")
    assert json.loads(result.stdout) == ["FutureWarning"]


def test_forwarding_modules_are_distinct_but_public_objects_are_identical():
    result = _fresh_python("""
import json
import bse.tools
import bse.point_group_data.C1
import chiq.tools
import chiq.point_group_data.C1
print(json.dumps({
    'tools_modules_distinct': bse.tools is not chiq.tools,
    'tools_object_identical': bse.tools.WallTime is chiq.tools.WallTime,
    'point_modules_distinct': bse.point_group_data.C1 is not chiq.point_group_data.C1,
    'point_object_identical': bse.point_group_data.C1.PointGroupData is chiq.point_group_data.C1.PointGroupData,
}))
""")
    assert json.loads(result.stdout) == {
        "tools_modules_distinct": True,
        "tools_object_identical": True,
        "point_modules_distinct": True,
        "point_object_identical": True,
    }


def test_all_allowlisted_modules_importable():
    for mod in SHIM_MODULES:
        assert importlib.util.find_spec("bse.%s" % mod) is not None
    for mod in ("h5bse", "bse_toml", "matrix_dict", "point_group", "tools", "index_pair", "mpi"):
        importlib.import_module("bse.%s" % mod)
    for name in POINT_GROUP_MODULES:
        importlib.import_module("bse.point_group_data.%s" % name)
