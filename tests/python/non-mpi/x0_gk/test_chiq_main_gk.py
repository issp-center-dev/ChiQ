"""chiq_main gives the same results from a G(k) file as from the equivalent saved-X0(q) file."""

import os
import shutil
import subprocess
import sys
from itertools import product

import h5py
import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(__file__))
import synth_gk as S  # noqa: E402

import chiq  # noqa: E402
from chiq import h5bse, x0_gk  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
TEMPLATE = os.path.join(HERE, "..", "bsetool_BSE", "ref", "dmft_bse.h5")
PACKAGE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(chiq.__file__)))

DIV = (4, 4, 1)
NUM_WF = 10  # matches the template's X_loc
TARGETS = ["00.00.00", "01.00.00", "02.01.00", "03.03.00", "01.02.00"]


def _run(module, args, cwd):
    env = dict(os.environ)
    env["PYTHONPATH"] = PACKAGE_DIR + os.pathsep + env.get("PYTHONPATH", "")
    return subprocess.run([sys.executable, "-m", module] + args, cwd=cwd, env=env,
                          capture_output=True, text=True)


def _template_layout():
    with h5py.File(TEMPLATE, "r") as f:
        beta = float(f["bse/info/beta"][()])
        names = [f["bse/info/block_name/" + str(i)][()].decode() for i in range(len(f["bse/info/block_name"]))]
    spins = ["up", "down"]
    block2 = []
    for name in names:  # e.g. '0-up-0-down' -> (c1, b1, c2, b2)
        c1, s1, c2, s2 = name.split("-")
        block2.append((int(c1), spins.index(s1), int(c2), spins.index(s2)))
    return beta, block2, [(0, 0)]


def _base_copy(path, phi):
    shutil.copy(TEMPLATE, path)
    with h5py.File(path, "a") as f:
        for name in ("X0_q", "X0_loc", "chi0_loc_in"):
            if "bse/input/" + name in f:
                del f["bse/input/" + name]
    h5bse.h5BSE(path).save(key=("Phi", 0), data=phi)


@pytest.fixture(scope="module")
def files(tmp_path_factory):
    work = tmp_path_factory.mktemp("gk_e2e")
    beta, block2, inner2 = _template_layout()
    n_iwG = 2 * NUM_WF
    G = S.model_g(beta, n_b1=2, n_c=1, n_in=1, n_iwG=n_iwG, div=DIV, seed=5)
    G = [0.3 * g for g in G]  # keep X0 small compared with the template's X_loc
    rng = np.random.default_rng(7)
    phi = {(b, b): 0.01 * (rng.normal(size=(1, 1, 2 * NUM_WF)) + 0j) for b in range(len(block2))}

    # saved file: X0_q for all q and X0_loc from the legacy FFT formula
    saved = str(work / "saved.h5")
    _base_copy(saved, phi)
    fft = S.fft_x0_all_q(G, block2, inner2, 0, NUM_WF, n_iwG // 2, DIV)
    bs = h5bse.h5BSE(saved)
    for iq, q in enumerate(np.ndindex(*DIV)):
        label = "%02d.%02d.%02d" % q
        bs.save(key=("X0_q", 0, label), data={k: v[..., iq] for k, v in fft.items()})
    bs.save(key=("X0_loc", 0), data={k: v.mean(axis=-1) for k, v in fft.items()})

    # gk file: x0_source and X0_loc from k-averaged G
    gk = str(work / "gk.h5")
    _base_copy(gk, phi)
    with h5py.File(gk, "a") as f:
        w = x0_gk.X0SourceWriter(f)
        w.write_meta(["up", "down"], block2, inner2, NUM_WF, [0], n_iwG, DIV, beta, (1, 1))
        for b in range(2):
            w.write_g_slice(b, 0, 0, 0, 0, G[b][0, 0, 0, 0])
        gbar = w.Gbar
        w.finish()
    h5bse.h5BSE(gk).save(key=("X0_loc", 0),
                         data=x0_gk.x0_loc_from_gbar(gbar, block2, inner2, 0, NUM_WF, n_iwG // 2))

    with open(work / "q.dat", "w") as f:
        for q in TARGETS:
            f.write("0 %s\n" % q)
    return work


def _toml(work, name, infile, outfile, types, omega_q="q.dat", extra=""):
    path = work / name
    path.write_text(
        "[chiq_common]\ninput = %r\noutput = %r\ntype = %r\n%s\n[chiq_main]\nbackend = \"numpy\"\n%s\n"
        % (infile, outfile, types, ("omega_q = %r" % omega_q) if omega_q else "", extra)
    )
    return str(path)


def _outputs(path):
    out = {}
    with h5py.File(path, "r") as f:
        def visit(name, obj):
            if isinstance(obj, h5py.Dataset) and name.startswith("bse/output/"):
                out[name] = obj[()]
        f.visititems(visit)
    return out


@pytest.mark.parametrize("num_wf", [None, 6])
def test_gk_and_saved_give_same_results(files, num_wf):
    work = files
    extra = "" if num_wf is None else "num_wf = %d" % num_wf
    types = ["chi0", "bse", "scl", "rpa", "rrpa"]
    results = {}
    for mode in ("saved", "gk"):
        out = "%s_%s.out.h5" % (mode, num_wf)
        toml = _toml(work, "%s_%s.toml" % (mode, num_wf), "%s.h5" % mode, out, types, extra=extra)
        r = _run("chiq.cli.chiq_main", [toml], cwd=str(work))
        assert r.returncode == 0, r.stdout + r.stderr
        results[mode] = _outputs(str(work / out))

    saved, gk = results["saved"], results["gk"]
    assert saved and set(saved) == set(gk)
    for kind in ("chi0_q", "chi_q", "chi_q_scl", "chi_q_rpa", "chi_q_rrpa"):
        assert any("/%s/" % kind in name for name in saved), kind
    # Tolerance relative to the scale of each quantity: I_q = X0^-1 - X^-1 has small
    # elements produced by cancellation, where elementwise relative error is not meaningful.
    scale = {}
    for name, value in saved.items():
        kind = name.split("/")[2]
        scale[kind] = max(scale.get(kind, 0.0), float(np.max(np.abs(value))))
    for name in saved:
        kind = name.split("/")[2]
        np.testing.assert_allclose(gk[name], saved[name], rtol=1e-10, atol=1e-10 * scale[kind], err_msg=name)


def test_gk_computes_exactly_the_requested_targets(tmp_path):
    rng = np.random.default_rng(3)
    div, num_wf, list_wb = (3, 2, 1), 2, [0, 1]
    n_iwG = 2 * (num_wf + 1)
    G = S.random_g(rng, 1, 1, 1, n_iwG, 6)
    block2, inner2 = S.make_pairs(1, 1, 1)
    S.write_gk_file(str(tmp_path / "in.h5"), G, block2, inner2, num_wf, list_wb, div, 5.0)
    (tmp_path / "q.dat").write_text("0 00.00.00\n1 02.01.00\n")
    toml = _toml(tmp_path, "t.toml", "in.h5", "out.h5", ["chi0"])
    r = _run("chiq.cli.chiq_main", [toml], cwd=str(tmp_path))
    assert r.returncode == 0, r.stdout + r.stderr
    with h5py.File(str(tmp_path / "out.h5"), "r") as f:
        got = {(w, q) for w in f["bse/output/chi0_q"] for q in f["bse/output/chi0_q/" + w]}
    assert got == {("w0", "q_00.00.00"), ("w1", "q_02.01.00")}


@pytest.mark.parametrize("line", ["0 03.00.00", "2 00.00.00"])
def test_gk_rejects_invalid_targets(tmp_path, line):
    rng = np.random.default_rng(4)
    G = S.random_g(rng, 1, 1, 1, 4, 6)
    block2, inner2 = S.make_pairs(1, 1, 1)
    S.write_gk_file(str(tmp_path / "in.h5"), G, block2, inner2, 2, [0], (3, 2, 1), 5.0)
    (tmp_path / "q.dat").write_text(line + "\n")
    toml = _toml(tmp_path, "t.toml", "in.h5", "out.h5", ["chi0"])
    r = _run("chiq.cli.chiq_main", [toml], cwd=str(tmp_path))
    assert r.returncode != 0
    assert "Error in omega_q" in r.stderr


def test_ambiguous_file_is_rejected(files, tmp_path):
    path = str(tmp_path / "both.h5")
    shutil.copy(str(files / "gk.h5"), path)
    with h5py.File(path, "a") as f:
        f.create_group("bse/input/X0_q")
    shutil.copy(str(files / "q.dat"), str(tmp_path / "q.dat"))
    toml = _toml(tmp_path, "t.toml", "both.h5", "out.h5", ["chi0"])
    r = _run("chiq.cli.chiq_main", [toml], cwd=str(tmp_path))
    assert r.returncode != 0
    assert "both 'X0_q' and 'x0_source'" in r.stderr


def test_calc_iq_scl_rejects_gk_file(files, tmp_path):
    shutil.copy(str(files / "gk.h5"), str(tmp_path / "dmft_bse.h5"))
    scl_in = os.path.join(HERE, "..", "bsetool_SCL", "scl_2pole.in")
    shutil.copy(scl_in, str(tmp_path / "scl_2pole.in"))
    r = _run("chiq.cli.calc_Iq_scl", ["scl_2pole.in"], cwd=str(tmp_path))
    assert r.returncode != 0
    assert "calc_Iq_scl does not support this yet" in r.stderr
