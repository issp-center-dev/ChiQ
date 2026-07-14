import json
import os
from pathlib import Path
import subprocess
import sys

import pytest


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


def test_verifier_executes_no_python_before_empty_environment_boundary():
    script = (ROOT / "tests/python/non-mpi/packaging/verify_pip_install.sh").read_text()
    wrapper = script.split("exec /usr/bin/env -i", 1)[0]
    assert " -c " not in wrapper
    assert "command -v" not in wrapper
    assert "dirname" not in wrapper
    assert "[[" not in wrapper
    assert "clean=(" not in wrapper


def test_verifier_normal_wrapper_scrubs_sitecustomize_before_python(tmp_path):
    script = ROOT / "tests/python/non-mpi/packaging/verify_pip_install.sh"
    poison = tmp_path / "poison"
    poison.mkdir()
    marker = tmp_path / "sitecustomize-ran"
    startup_marker = tmp_path / "bash-env-ran"
    bash_env = tmp_path / "bash-env.sh"
    bash_env.write_text("/usr/bin/touch %s\n" % startup_marker)
    (poison / "sitecustomize.py").write_text(
        "from pathlib import Path\nPath(%r).write_text('ran')\n" % str(marker)
    )
    env = {
        "HOME": os.environ.get("HOME", str(tmp_path)),
        "PATH": str(poison),
        "PYTHON": sys.executable,
        "PYTHONPATH": str(poison),
        "PIP_INDEX_URL": "https://poison.invalid/simple",
        "GIT_DIR": str(tmp_path / "poison.git"),
        "BASH_ENV": str(bash_env),
        "TMPDIR": str(tmp_path),
    }

    result = subprocess.run(
        [str(script), "--audit-environment-only"],
        cwd=str(ROOT),
        env=env,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, result.stderr
    assert "Environment policy: env -i explicit allowlist" in result.stdout
    assert "Environment audit OK" in result.stdout
    assert not marker.exists()
    assert not startup_marker.exists()


@pytest.mark.parametrize("through_symlink", [False, True])
def test_verifier_rehomes_tmpdir_outside_source_tree(tmp_path, through_symlink):
    script = ROOT / "tests/python/non-mpi/packaging/verify_pip_install.sh"
    if through_symlink:
        link = tmp_path / "source-link"
        link.symlink_to(ROOT, target_is_directory=True)
        requested = link / "nested-diagnostics"
    else:
        requested = ROOT / "nested-diagnostics"
    env = {
        "HOME": os.environ.get("HOME", str(tmp_path)),
        "PATH": "/usr/bin:/bin",
        "PYTHON": sys.executable,
        "TMPDIR": str(requested),
    }

    result = subprocess.run(
        [str(script), "--audit-environment-only"],
        cwd=str(ROOT),
        env=env,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, result.stderr
    assert "TMPDIR=%s\n" % Path("/tmp").resolve() in result.stdout
    assert str(ROOT) not in next(
        line for line in result.stdout.splitlines() if line.startswith("TMPDIR=")
    )


def test_verifier_rejects_poisoned_direct_internal_mode_before_work(tmp_path):
    script = ROOT / "tests/python/non-mpi/packaging/verify_pip_install.sh"
    poison = tmp_path / "poison"
    poison.mkdir()
    env = {
        "CHIQ_VERIFIER_INTERNAL": "1",
        "PATH": str(poison),
        "PYTHON_REQUEST": sys.executable,
        "PYTHONPATH": str(poison),
        "PIP_CONFIG_FILE": str(tmp_path / "pip.conf"),
        "GIT_WORK_TREE": str(tmp_path / "checkout"),
        "BASH_ENV": str(tmp_path / "missing-startup-file"),
    }

    result = subprocess.run(
        ["/bin/bash", str(script), "--audit-environment-only"],
        cwd=str(ROOT),
        env=env,
        capture_output=True,
        text=True,
    )

    assert result.returncode != 0
    assert "unsafe internal verifier environment" in result.stderr
    assert "Environment policy" not in result.stdout
    assert "source snapshot" not in result.stdout + result.stderr
