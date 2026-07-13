import importlib
import os
import subprocess
import sys
from types import SimpleNamespace

import pytest

from contracts import CANONICAL_COMMANDS


ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(__file__)))))


def _core_only_env(tmp_path):
    blocker = tmp_path / "blocked_optional_dependencies"
    blocker.mkdir()
    (blocker / "matplotlib.py").write_text("raise ImportError('matplotlib blocked')\n")
    (blocker / "dcore.py").write_text("raise ImportError('dcore blocked')\n")
    env = dict(os.environ)
    env["PYTHONPATH"] = os.pathsep.join((
        str(blocker),
        os.path.join(ROOT, "python", "package"),
        os.path.join(ROOT, "build", "src"),
    ))
    return env


def _run_module(name, args, tmp_path):
    code = "from chiq.cli.%s import main; main()" % name
    return subprocess.run(
        [sys.executable, "-c", code] + list(args),
        cwd=str(tmp_path),
        env=_core_only_env(tmp_path),
        capture_output=True,
        text=True,
    )


@pytest.mark.parametrize("name", CANONICAL_COMMANDS)
def test_cli_module_has_main(name):
    mod = importlib.import_module("chiq.cli.%s" % name)
    assert callable(mod.main)

def test_gen_qpath_exports_GenQPath():
    from chiq.cli.gen_qpath import GenQPath
    assert GenQPath is not None

@pytest.mark.parametrize("name", CANONICAL_COMMANDS)
def test_canonical_help_runs_on_core_install(name, tmp_path):
    r = _run_module(name, ["--help"], tmp_path)
    assert r.returncode == 0, r.stderr
    assert "usage:" in r.stdout.lower()
    assert "Traceback" not in r.stderr


@pytest.mark.parametrize("name", CANONICAL_COMMANDS)
def test_canonical_version_is_exact_and_runs_on_core_install(name, tmp_path):
    from chiq import __version__

    r = _run_module(name, ["--version"], tmp_path)
    assert r.returncode == 0, r.stderr
    assert r.stdout == "ChiQ version %s\n" % __version__
    assert r.stderr == ""
    assert "Traceback" not in r.stdout + r.stderr


@pytest.mark.parametrize(
    "name,args",
    (
        ("plot_chiq_path", ("q-path.dat", "eigen.dat")),
        ("plot_Ir", ("I-r.dat",)),
    ),
)
def test_plot_commands_report_missing_plot_extra_after_valid_parse(name, args, tmp_path):
    for filename in args:
        (tmp_path / filename).write_text("")
    r = _run_module(name, args, tmp_path)
    assert r.returncode != 0
    assert "pip install chiq[plot]" in r.stderr
    assert "Traceback" not in r.stdout + r.stderr


@pytest.mark.parametrize("name", CANONICAL_COMMANDS)
def test_scripts_are_thin_wrappers(name):
    p = os.path.join(ROOT, "python", "scripts", name + ".py")
    text = open(p).read()
    assert "from chiq.cli.%s import main" % name in text


@pytest.mark.parametrize("name", CANONICAL_COMMANDS)
def test_deprecated_alias_prints_stderr_and_forwards_argv(name, monkeypatch, capsys):
    from chiq.cli import _deprecated

    argv = [name + ".py", "--sentinel", "value"]
    seen = []

    def forwarded_main():
        seen.append(list(sys.argv))
        return "forwarded-result"

    monkeypatch.setattr(sys, "argv", argv)
    monkeypatch.setattr(
        _deprecated.importlib,
        "import_module",
        lambda module_name: SimpleNamespace(main=forwarded_main),
    )
    result = getattr(_deprecated, name + "_py")()

    captured = capsys.readouterr()
    assert captured.out == ""
    assert "'%s.py' command is deprecated" % name in captured.err
    assert "use '%s' instead" % name in captured.err
    assert seen == [argv]
    assert result == "forwarded-result"
