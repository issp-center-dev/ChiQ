import importlib.util
import os
from pathlib import Path
import subprocess
import sys

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


def test_boundary_smoke_is_python38_compatible():
    source = SMOKE_PATH.read_text()
    assert "list[" not in source
    assert "dict[" not in source
    assert " | None" not in source
    compile(source, str(SMOKE_PATH), "exec")


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
