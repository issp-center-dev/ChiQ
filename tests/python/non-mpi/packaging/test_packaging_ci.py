import json
from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[4]


def test_ci_runs_packaging_verification_scripts():
    workflow = (ROOT / ".github/workflows/main.yml").read_text()
    assert "verify_legacy_install.sh" in workflow
    assert "verify_pip_install.sh" in workflow
    assert "submodules: recursive" in workflow


def test_packaging_verification_scripts_support_linux():
    packaging = ROOT / "tests/python/non-mpi/packaging"
    for name in ("verify_legacy_install.sh", "verify_pip_install.sh"):
        script = (packaging / name).read_text()
        assert 'uname -s' in script
        assert "/usr/include/eigen3" in script


def test_pip_verifier_uses_safe_snapshot_and_artifact_chain():
    script = (ROOT / "tests/python/non-mpi/packaging/verify_pip_install.sh").read_text()
    assert "source_snapshot.py" in script
    assert "snapshot-manifest.json" in script
    assert "cp -R" not in script
    assert "artifact_inspector.py\" sdist" in script
    assert "--extract" in script
    assert "artifact_inspector.py\" wheel" in script
    assert "WHEEL_SOURCE=\"$EXTRACTED_ROOT\"" in script
    assert "-m build --wheel" in script
    assert "cd \"$WHEEL_SOURCE\"" in script
    assert "find \"$EXTRACTED_ROOT\" -name .git" in script
    assert "extern/pybind11" in script


def test_pip_verifier_has_four_distinct_venvs_and_exact_installs():
    script = (ROOT / "tests/python/non-mpi/packaging/verify_pip_install.sh").read_text()
    for name in ("FRONTEND_VENV", "WHEEL_VENV", "SDIST_VENV", "EDITABLE_VENV"):
        assert name + "=" in script
    assert script.count("-m venv") == 4
    assert 'pip install --no-cache-dir "$WHEEL[mpi]"' in script
    assert 'pip install --no-cache-dir "$SDIST[mpi]"' in script
    assert 'pip install --no-cache-dir -e "$SNAPSHOT[mpi]"' in script
    assert script.count("verify_hash ") == 2
    assert 'runtime_smoke.py" --mode wheel' in script
    assert 'runtime_smoke.py" --mode sdist' in script
    assert 'runtime_smoke.py" --mode editable' in script
    assert 'cd "$EXTERNAL_CWD"' in script


def test_pip_verifier_sanitizes_and_separates_mode_state():
    script = (ROOT / "tests/python/non-mpi/packaging/verify_pip_install.sh").read_text()
    assert "env -i" in script
    for forbidden in (
        "PYTHONPATH", "PYTHONHOME", "VIRTUAL_ENV", "CONDA_PREFIX",
        "CMAKE_PREFIX_PATH", "CMAKE_BUILD_PARALLEL_LEVEL", "SKBUILD_",
        "PIP_", "GIT_",
    ):
        assert forbidden in script
    for mode in ("wheel", "sdist", "editable"):
        assert 'TMPDIR="$DIAGNOSTICS/tmp/%s"' % mode in script
        assert 'PIP_CACHE_DIR="$DIAGNOSTICS/cache/%s"' % mode in script
    assert "sanitized-environment.txt" in script
    assert "frontend-versions.txt" in script
    assert "backend-distributions.json" in script
    assert "pip-list.json" in script
    assert "trap" in script and "Diagnostics retained at" in script


def test_pip_verifier_pins_python38_frontend():
    script = (ROOT / "tests/python/non-mpi/packaging/verify_pip_install.sh").read_text()
    assert '"pip<25.1"' in script
    assert '"build<1.3"' in script


def test_runtime_smokes_use_only_their_install_venv_commands():
    script = (ROOT / "tests/python/non-mpi/packaging/verify_pip_install.sh").read_text()
    for venv in ("WHEEL_VENV", "SDIST_VENV", "EDITABLE_VENV"):
        prefix = 'PATH="$%s/bin:$PATH"' % venv
        assert script.count(prefix) == 1
    assert 'PATH="$REPO' not in script
    assert 'PATH="$SNAPSHOT' not in script
    assert 'PATH="$REPO/build' not in script


def test_backend_probe_canonicalizes_pip_distribution_names(tmp_path):
    script = (ROOT / "tests/python/non-mpi/packaging/verify_pip_install.sh").read_text()
    parser = script.split("# BACKEND_PROBE_PARSER_BEGIN\n", 1)[1].split(
        "# BACKEND_PROBE_PARSER_END", 1
    )[0]
    source = tmp_path / "pip-list.json"
    destination = tmp_path / "backend-distributions.json"
    source.write_text(json.dumps([
        {"name": "scikit_build_core", "version": "0.12.2"},
        {"name": "pybind11", "version": "2.13.6"},
        {"name": "cmake", "version": "4.4.0"},
        {"name": "unrelated", "version": "1.0"},
    ]))

    result = subprocess.run(
        [sys.executable, "-c", parser, str(source), str(destination)],
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, result.stderr
    assert {row["name"] for row in json.loads(destination.read_text())} == {
        "scikit_build_core", "pybind11", "cmake"
    }
