import ast
import importlib.util
import os
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace

import pytest


ROOT = Path(__file__).resolve().parents[4]
SMOKE_PATH = Path(__file__).with_name("dcore_boundary_smoke.py")


def _core_only_environment(tmp_path):
    blocker = tmp_path / "blocked"
    blocker.mkdir()
    (blocker / "dcore.py").write_text("raise ImportError('dcore is not installed')\n")
    environment = dict(os.environ)
    environment["PYTHONPATH"] = os.pathsep.join((
        str(blocker),
        str(ROOT / "python" / "package"),
        str(ROOT / "build" / "src"),
    ))
    return environment


@pytest.mark.parametrize("deprecated", [False, True])
def test_missing_dcore_extra_is_actionable_without_traceback(tmp_path, deprecated):
    command_name = "dcore_chiq.py" if deprecated else "dcore_chiq"
    if deprecated:
        code = (
            "import sys; sys.argv[0] = %r; "
            "from chiq.cli._deprecated import dcore_chiq_py; dcore_chiq_py()"
            % command_name
        )
    else:
        code = (
            "import sys; sys.argv[0] = %r; "
            "from chiq.cli.dcore_chiq import main; main()" % command_name
        )
    result = subprocess.run(
        [sys.executable, "-c", code, str(tmp_path / "missing.ini"), "--np", "1"],
        cwd=str(tmp_path),
        env=_core_only_environment(tmp_path),
        capture_output=True,
        text=True,
    )

    assert result.returncode != 0
    assert "pip install chiq[dcore]" in result.stderr
    assert "Traceback" not in result.stdout + result.stderr
    if deprecated:
        assert "deprecated" in result.stderr.lower()


def test_boundary_smoke_checks_every_private_dcore_symbol_used_by_cli():
    spec = importlib.util.spec_from_file_location("dcore_boundary_smoke", str(SMOKE_PATH))
    assert spec is not None and spec.loader is not None
    smoke = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(smoke)

    assert smoke.REQUIRED_CALLABLES == {
        "dcore._dispatcher": ("HDFArchive", "dyson"),
        "dcore.dmft_core": ("DMFTCoreSolver",),
        "dcore.program_options": (
            "create_parser", "parse_parameters", "print_parameters",
            "delete_parameters", "OptionStatus",
        ),
        "dcore.tools": (
            "raise_if_mpi_imported", "make_empty_dir",
            "launch_mpi_subprocesses", "float_to_complex_array",
        ),
        "dcore.sumkdft_workers.launcher": ("run_sumkdft",),
    }
    assert smoke.REQUIRED_MEMBERS == {
        "dcore.impurity_solvers": ("solver_classes", "compute_basis_rot"),
        "mpi4py": ("MPI",),
    }


def _load_smoke_module():
    spec = importlib.util.spec_from_file_location("dcore_boundary_smoke", str(SMOKE_PATH))
    assert spec is not None and spec.loader is not None
    smoke = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(smoke)
    return smoke


def _fake_boundary_modules(root):
    modules = {
        "dcore": SimpleNamespace(__file__=str(root / "dcore" / "__init__.py")),
    }
    for module_name, names in _load_smoke_module().REQUIRED_CALLABLES.items():
        values = {name: (lambda: None) for name in names}
        values["__file__"] = str(root / (module_name.replace(".", "/") + ".py"))
        modules[module_name] = SimpleNamespace(**values)
    modules["dcore.impurity_solvers"] = SimpleNamespace(
        __file__=str(root / "dcore" / "impurity_solvers.py"),
        solver_classes={},
        compute_basis_rot=lambda: None,
    )
    modules["mpi4py"] = SimpleNamespace(
        __file__=str(root / "mpi4py" / "__init__.py"),
        MPI=object(),
    )
    return modules


def test_installed_boundary_checks_dependency_versions_and_all_module_origins(
        tmp_path, monkeypatch):
    smoke = _load_smoke_module()
    expected_root = tmp_path / "venv"
    modules = _fake_boundary_modules(expected_root)
    versions = {"dcore": "4.2.0", "mpi4py": "4.1.1"}
    seen_versions = []

    monkeypatch.setattr(smoke.importlib, "import_module", modules.__getitem__)

    def version(name):
        seen_versions.append(name)
        return versions[name]

    monkeypatch.setattr(smoke.importlib.metadata, "version", version)
    smoke._check_import_boundary(expected_root)

    assert seen_versions == ["dcore", "mpi4py"]


@pytest.mark.parametrize("outside", ["mpi4py", "dcore.tools"])
def test_installed_boundary_rejects_dependency_module_outside_environment(
        tmp_path, monkeypatch, outside):
    smoke = _load_smoke_module()
    expected_root = tmp_path / "venv"
    modules = _fake_boundary_modules(expected_root)
    modules[outside].__file__ = str(tmp_path / "foreign" / "module.py")
    monkeypatch.setattr(smoke.importlib, "import_module", modules.__getitem__)
    monkeypatch.setattr(
        smoke.importlib.metadata,
        "version",
        lambda name: {"dcore": "4.2.0", "mpi4py": "4.1.1"}[name],
    )

    with pytest.raises(smoke.BoundaryError, match="outside the install environment"):
        smoke._check_import_boundary(expected_root)


def test_installed_boundary_reports_module_without_filesystem_origin(
        tmp_path, monkeypatch):
    smoke = _load_smoke_module()
    expected_root = tmp_path / "venv"
    modules = _fake_boundary_modules(expected_root)
    del modules["mpi4py"].__file__
    monkeypatch.setattr(smoke.importlib, "import_module", modules.__getitem__)
    monkeypatch.setattr(
        smoke.importlib.metadata,
        "version",
        lambda name: {"dcore": "4.2.0", "mpi4py": "4.1.1"}[name],
    )

    with pytest.raises(smoke.BoundaryError, match="has no filesystem origin"):
        smoke._check_import_boundary(expected_root)


def test_boundary_smoke_is_python38_compatible():
    source = SMOKE_PATH.read_text()
    assert "list[" not in source
    assert "dict[" not in source
    assert " | None" not in source
    try:
        tree = ast.parse(source, filename=str(SMOKE_PATH), feature_version=(3, 8))
    except TypeError:
        tree = ast.parse(source, filename=str(SMOKE_PATH), feature_version=8)
    assert not any(isinstance(node, ast.Assert) for node in ast.walk(tree))


def test_optimized_smoke_rejects_foreign_origin(tmp_path):
    code = (
        "import importlib.util, sys; "
        "spec = importlib.util.spec_from_file_location('smoke', sys.argv[1]); "
        "smoke = importlib.util.module_from_spec(spec); "
        "spec.loader.exec_module(smoke); "
        "smoke._assert_installed_origin('mpi4py', sys.argv[2], sys.argv[3])"
    )
    result = subprocess.run(
        [
            sys.executable, "-O", "-c", code, str(SMOKE_PATH),
            str(tmp_path / "foreign" / "mpi4py.py"), str(tmp_path / "venv"),
        ],
        capture_output=True,
        text=True,
    )

    assert result.returncode != 0
    assert "outside the install environment" in result.stderr


def test_smoke_main_reports_controlled_boundary_failure(monkeypatch, capsys):
    smoke = _load_smoke_module()

    def fail(_argv):
        raise smoke.BoundaryError("sentinel boundary failure")

    monkeypatch.setattr(smoke, "_main", fail)
    assert smoke.main([]) == 1
    captured = capsys.readouterr()
    assert "sentinel boundary failure" in captured.err
    assert "Traceback" not in captured.err


def test_ci_has_endpoint_dcore_failure_domains_and_retains_diagnostics():
    workflow = (ROOT / ".github/workflows/main.yml").read_text()
    for job in ("packaging-unit:", "pip-artifacts:", "core-boundary:", "dcore-boundary:"):
        assert job in workflow
    assert workflow.count("python-version: ['3.8', '3.13']") >= 3
    assert "openmpi-bin libopenmpi-dev" in workflow
    assert "pip install --no-cache-dir '.[dcore]'" in workflow
    assert "pip install --no-cache-dir ." in workflow
    assert "PIP_NO_CACHE_DIR: \"1\"" in workflow
    assert 'dcore_boundary_smoke.py" --mode installed' in workflow
    assert 'dcore_boundary_smoke.py" --mode missing-extra' in workflow
    assert "actions/upload-artifact@v4" in workflow
    assert "if: failure()" in workflow
    assert "chiq-pip-verification.*" in workflow


def test_ci_boundary_install_and_smoke_pipelines_fail_closed_and_retain_logs():
    workflow = (ROOT / ".github/workflows/main.yml").read_text()
    core = workflow.split("  core-boundary:", 1)[1].split(
        "\n  dcore-boundary:", 1
    )[0]
    dcore = workflow.split("  dcore-boundary:", 1)[1].split("\n  legacy:", 1)[0]

    for block, mode, log_prefix in (
        (core, "missing-extra", "core-boundary"),
        (dcore, "installed", "dcore-boundary"),
    ):
        assert block.count("set -o pipefail") >= 2
        smoke_pipeline = next(
            line for line in block.splitlines()
            if "dcore_boundary_smoke.py" in line
        )
        assert "--mode %s" % mode in smoke_pipeline
        assert '| tee "$RUNNER_TEMP/%s.log"' % log_prefix in smoke_pipeline
        assert 'tee "$RUNNER_TEMP/%s-install.log"' % log_prefix in block
        assert "${{ runner.temp }}/%s-install.log" % log_prefix in block
    assert "if-no-files-found: warn" in dcore


def test_ci_retains_endpoint_scientific_and_cpp_regressions():
    workflow = (ROOT / ".github/workflows/main.yml").read_text()
    scientific = workflow.split("  scientific:", 1)[1].split(
        "\n  packaging-unit:", 1
    )[0]
    assert "python-version: ['3.8', '3.13']" in scientific
    assert "cmake .. -DTesting=ON" in scientific
    assert "ctest -V" in scientific
    assert "python -m pytest -v" in scientific


def test_ci_runs_authoritative_pip_flow_once_per_endpoint():
    workflow = (ROOT / ".github/workflows/main.yml").read_text()
    pip_job = workflow.split("  pip-artifacts:", 1)[1].split("\n  core-boundary:", 1)[0]
    assert "python-version: ['3.8', '3.13']" in pip_job
    assert pip_job.count("verify_pip_install.sh") == 1
    assert "timeout-minutes:" in pip_job


def test_authoritative_flow_uses_endpoint_build_backend_requirement():
    script = (
        ROOT / "tests/python/non-mpi/packaging/verify_pip_install.sh"
    ).read_text()
    assert 'BACKEND_REQUIREMENT="scikit-build-core>=0.10,<0.11"' in script
    assert 'BACKEND_REQUIREMENT="scikit-build-core>=0.10,<1"' in script
    assert '"$BACKEND_REQUIREMENT"' in script
