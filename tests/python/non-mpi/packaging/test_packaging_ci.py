from pathlib import Path


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
