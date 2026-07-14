import os
import sys

from contracts import CANONICAL_COMMANDS, DEPRECATED_COMMANDS


if sys.version_info >= (3, 11):
    import tomllib

    def _load(path):
        with open(path, "rb") as f:
            return tomllib.load(f)
else:
    import toml

    def _load(path):
        return toml.load(path)


ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(__file__)))))


def test_pyproject_core_metadata():
    d = _load(os.path.join(ROOT, "pyproject.toml"))
    assert d["build-system"]["build-backend"] == "scikit_build_core.build"
    assert d["project"]["requires-python"] == ">=3.8"
    core = " ".join(d["project"]["dependencies"])
    for dep in ["numpy", "scipy", "more-itertools", "h5py", "toml"]:
        assert dep in core
    extras = d["project"]["optional-dependencies"]
    assert "matplotlib" in " ".join(extras["plot"])
    assert "mpi4py" in " ".join(extras["mpi"])
    assert "dcore" in " ".join(extras["dcore"])


def test_pyproject_pins_endpoint_build_frontends_and_dcore_boundary():
    d = _load(os.path.join(ROOT, "pyproject.toml"))
    assert d["build-system"]["requires"] == [
        "scikit-build-core>=0.10,<0.11; python_version < '3.9'",
        "scikit-build-core>=0.10,<1; python_version >= '3.9'",
        "pybind11>=2.12,<3",
        "cmake>=3.15",
    ]
    extras = d["project"]["optional-dependencies"]
    assert extras["dcore"] == ["dcore==4.2.0", "mpi4py"]
    assert extras["mpi"] == ["mpi4py"]


def test_pyproject_has_exact_script_contract():
    d = _load(os.path.join(ROOT, "pyproject.toml"))
    scripts = d["project"]["scripts"]
    expected = {}
    for name in CANONICAL_COMMANDS:
        expected[name] = "chiq.cli.%s:main" % name
        expected[name + ".py"] = "chiq.cli._deprecated:%s_py" % name
    assert len(scripts) == 22
    assert tuple(sorted(scripts)) == tuple(sorted(CANONICAL_COMMANDS + DEPRECATED_COMMANDS))
    assert scripts == expected


def test_pyproject_retains_packages_and_dynamic_version_provider():
    d = _load(os.path.join(ROOT, "pyproject.toml"))
    assert d["project"]["dynamic"] == ["version"]
    assert d["tool"]["scikit-build"]["wheel"]["packages"] == [
        "python/package/chiq",
        "python/package/bse",
        "python/package/bse_solver",
    ]
    version = d["tool"]["scikit-build"]["metadata"]["version"]
    assert version["provider"] == "scikit_build_core.metadata.regex"
    assert version["input"] == "python/package/chiq/__init__.py"
    assert "__version__" in version["regex"]


def test_pyproject_contains_no_local_macos_paths():
    text = open(os.path.join(ROOT, "pyproject.toml")).read()
    forbidden = (
        "/Users/",
        "/Library/Developer/",
        "/Applications/",
        "/opt/homebrew/",
        "/usr/local/opt/",
        "MacOSX.sdk",
    )
    assert not any(path in text for path in forbidden)


def test_pyproject_excludes_source_artifacts_without_hiding_built_extension():
    d = _load(os.path.join(ROOT, "pyproject.toml"))
    config = d["tool"]["scikit-build"]
    assert config["wheel"]["packages"] == [
        "python/package/chiq",
        "python/package/bse",
        "python/package/bse_solver",
    ]
    sdist_exclude = " ".join(config["sdist"]["exclude"])
    for pattern in ("_bse_solver", "__pycache__", "*.pyc", "build"):
        assert pattern in sdist_exclude
    assert not any("tests" in pattern and "*.h5" in pattern for pattern in config["sdist"]["exclude"])
    wheel_exclude = " ".join(config["wheel"]["exclude"])
    assert "__pycache__" in wheel_exclude
    assert "*.pyc" in wheel_exclude
    assert "_bse_solver" not in wheel_exclude
