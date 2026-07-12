import os
import sys


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
    scripts = d["project"]["scripts"]
    for cmd in ["chiq_main", "chiq_post", "gen_qpath", "dcore_chiq"]:
        assert scripts[cmd] == f"chiq.cli.{cmd}:main"
    assert scripts["chiq_main.py"] == "chiq.cli._deprecated:chiq_main_py"
    extras = d["project"]["optional-dependencies"]
    assert "matplotlib" in " ".join(extras["plot"])
    assert "mpi4py" in " ".join(extras["mpi"])
    assert "dcore" in " ".join(extras["dcore"])
