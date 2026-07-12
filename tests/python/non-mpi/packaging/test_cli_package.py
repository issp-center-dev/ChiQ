import importlib
import subprocess
import sys
import pytest

CORE_COMMANDS = [
    "chiq_main", "chiq_post", "chiq_fft", "gen_qpath", "gen_allq",
    "calc_Iq", "calc_Iq_scl", "plot_chiq_path", "plot_Ir", "eigenvec_viewer",
]

@pytest.mark.parametrize("name", CORE_COMMANDS + ["dcore_chiq"])
def test_cli_module_has_main(name):
    mod = importlib.import_module(f"chiq.cli.{name}")
    assert callable(mod.main)

def test_gen_qpath_exports_GenQPath():
    from chiq.cli.gen_qpath import GenQPath
    assert GenQPath is not None

def test_chiq_main_version_runs_on_core():
    # --version must work without optional deps / without the built extension.
    # `python -c "<code>" --version` puts "--version" in sys.argv[1], which argparse reads.
    import os
    env = dict(os.environ)
    r = subprocess.run([sys.executable, "-c",
                        "from chiq.cli.chiq_main import main; main()", "--version"],
                       env=env, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    assert "ChiQ version" in (r.stdout + r.stderr)
