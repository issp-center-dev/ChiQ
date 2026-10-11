import os
import sys

import h5py
import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(__file__))
import synth_gk as S  # noqa: E402

from chiq import x0_gk  # noqa: E402
from chiq.x0_gk import X0SourceError, GkX0Source  # noqa: E402

DIV = (4, 3, 2)
N_K = 24
NUM_WF = 2
LIST_WB = [0, 1]
N_IWG = 2 * (NUM_WF + max(LIST_WB))
BETA = 7.5


def _assert_dict_close(a, b, rtol=1e-12, atol=1e-14):
    assert set(a) == set(b)
    for key in a:
        np.testing.assert_allclose(a[key], b[key], rtol=rtol, atol=atol)


@pytest.fixture(scope="module")
def setup():
    rng = np.random.default_rng(1)
    n_b1, n_c, n_in = 2, 2, 2
    G = S.random_g(rng, n_b1, n_c, n_in, N_IWG, N_K)
    block2, inner2 = S.make_pairs(n_c, n_b1, n_in, nonlocal_pairs=True)
    return G, block2, inner2


@pytest.mark.parametrize("wb", LIST_WB)
@pytest.mark.parametrize("q", [(0, 0, 0), (1, 0, 0), (3, 2, 1)])
def test_kernel_matches_loop_oracle(setup, wb, q):
    G, block2, inner2 = setup
    got = x0_gk.x0_from_g(G, block2, inner2, wb, NUM_WF, N_IWG // 2, x0_gk.k_minus_q_index(DIV, q))
    ref = S.oracle_x0(G, block2, inner2, wb, NUM_WF, N_IWG // 2, DIV, q, sign=-1)
    _assert_dict_close(got, ref)


def test_sign_is_k_minus_q(setup):
    G, block2, inner2 = setup
    q = (1, 0, 0)
    got = x0_gk.x0_from_g(G, block2, inner2, 1, NUM_WF, N_IWG // 2, x0_gk.k_minus_q_index(DIV, q))
    plus = S.oracle_x0(G, block2, inner2, 1, NUM_WF, N_IWG // 2, DIV, q, sign=+1)
    key = next(iter(got))
    assert not np.allclose(got[key], plus[key])


@pytest.mark.parametrize("wb", LIST_WB)
def test_kernel_matches_legacy_fft_formula_for_all_q(setup, wb):
    G, block2, inner2 = setup
    ref = S.fft_x0_all_q(G, block2, inner2, wb, NUM_WF, N_IWG // 2, DIV)
    for iq, q in enumerate(np.ndindex(*DIV)):
        got = x0_gk.x0_from_g(G, block2, inner2, wb, NUM_WF, N_IWG // 2, x0_gk.k_minus_q_index(DIV, q))
        for key in got:
            np.testing.assert_allclose(got[key], ref[key][..., iq], rtol=1e-10, atol=1e-12)


@pytest.mark.parametrize("wb", LIST_WB)
def test_x0_loc_equals_average_over_all_q(setup, wb):
    G, block2, inner2 = setup
    gbar = [g.mean(axis=-1) for g in G]
    got = x0_gk.x0_loc_from_gbar(gbar, block2, inner2, wb, NUM_WF, N_IWG // 2)
    ref = S.fft_x0_all_q(G, block2, inner2, wb, NUM_WF, N_IWG // 2, DIV)
    for key in got:
        np.testing.assert_allclose(got[key], ref[key].mean(axis=-1), rtol=1e-10, atol=1e-12)


def test_soc_like_single_block():
    rng = np.random.default_rng(2)
    G = S.random_g(rng, 1, 1, 4, N_IWG, N_K)
    block2, inner2 = S.make_pairs(1, 1, 4)
    q = (2, 1, 1)
    got = x0_gk.x0_from_g(G, block2, inner2, 1, NUM_WF, N_IWG // 2, x0_gk.k_minus_q_index(DIV, q))
    ref = S.oracle_x0(G, block2, inner2, 1, NUM_WF, N_IWG // 2, DIV, q)
    _assert_dict_close(got, ref)


def test_local_block_pairs_only():
    rng = np.random.default_rng(3)
    G = S.random_g(rng, 2, 2, 1, N_IWG, N_K)
    block2, inner2 = S.make_pairs(2, 2, 1, nonlocal_pairs=False)
    q = (1, 2, 0)
    got = x0_gk.x0_from_g(G, block2, inner2, 0, NUM_WF, N_IWG // 2, x0_gk.k_minus_q_index(DIV, q))
    ref = S.oracle_x0(G, block2, inner2, 0, NUM_WF, N_IWG // 2, DIV, q)
    _assert_dict_close(got, ref)


# ---------------------------------------------------------------- provider

@pytest.fixture
def gk_file(setup, tmp_path):
    G, block2, inner2 = setup
    path = str(tmp_path / "gk.h5")
    S.write_gk_file(path, G, block2, inner2, NUM_WF, LIST_WB, DIV, BETA)
    return path


def test_provider_round_trip(setup, gk_file):
    G, block2, inner2 = setup
    src = GkX0Source(gk_file)
    for wb in LIST_WB:
        got = src.get(wb, "03.02.01")
        ref = S.oracle_x0(G, block2, inner2, wb, NUM_WF, N_IWG // 2, DIV, (3, 2, 1))
        _assert_dict_close(got, ref)
    with h5py.File(gk_file, "r") as h5f:
        assert x0_gk.detect_representation(h5f) == "gk"


@pytest.mark.parametrize("label", ["04.00.00", "00.03.00", "-1.00.00", "00.00", "a.b.c"])
def test_provider_rejects_q_off_mesh(gk_file, label):
    with pytest.raises(X0SourceError):
        GkX0Source(gk_file).get(0, label)


def test_provider_rejects_unknown_wb(gk_file):
    with pytest.raises(X0SourceError):
        GkX0Source(gk_file).get(2, "00.00.00")


def _mutate(path, func):
    with h5py.File(path, "a") as h5f:
        func(h5f["bse/x0_source"], h5f)


def _replace(grp, name, value):
    del grp[name]
    grp.create_dataset(name, data=value)


@pytest.mark.parametrize("mutation", [
    lambda g, f: g.attrs.__setitem__("complete", 0),
    lambda g, f: g.attrs.__setitem__("version", 99),
    lambda g, f: g.attrs.__setitem__("format", "other"),
    lambda g, f: f.create_group("bse/input/X0_q"),
    lambda g, f: g.attrs.__setitem__("beta", BETA * 2),
    lambda g, f: g.attrs.__setitem__("iw_offset", 1),
    lambda g, f: g.attrs.__setitem__("list_wb", np.array([0, 1, 5])),
    lambda g, f: _replace(g, "inner2_pairs", np.array([[0, 0], [0, -1], [1, 0], [1, 1]])),
    lambda g, f: _replace(g, "inner2_pairs", np.array([[0, 0], [0, 2], [1, 0], [1, 1]])),
    lambda g, f: _replace(g, "inner2_pairs", np.array([[0.0, 0.0], [0, 1], [1, 0], [1, 1]])),
    lambda g, f: _replace(g, "inner2_pairs", np.array([[0, 0, 0]])),
    lambda g, f: _replace(g, "block2_pairs", g["block2_pairs"][()][:-1]),
    lambda g, f: _replace(g, "block2_pairs", np.where(g["block2_pairs"][()] == 1, -1, g["block2_pairs"][()])),
    lambda g, f: _replace(g, "block2_pairs", np.where(g["block2_pairs"][()] == 1, 2, g["block2_pairs"][()])),
    lambda g, f: _replace(g, "G/1", g["G/1"][()][..., :-1]),
    lambda g, f: g.__delitem__("G/1"),
    lambda g, f: f.__delitem__("bse/input/X0_loc/w1"),
])
def test_provider_rejects_malformed_file(gk_file, mutation):
    _mutate(gk_file, mutation)
    with pytest.raises(X0SourceError):
        GkX0Source(gk_file)


def test_writer_replaces_saved_x0q(tmp_path):
    path = str(tmp_path / "f.h5")
    with h5py.File(path, "w") as h5f:
        h5f.create_dataset("bse/input/X0_q/w0/q_00.00.00/0_0", data=np.zeros((1, 1, 2)))
        h5f.create_dataset("bse/input/X_loc/w0/0_0", data=np.zeros((1, 1, 2, 2)))
        x0_gk.X0SourceWriter(h5f)
        assert "bse/input/X0_q" not in h5f
        assert "bse/input/X_loc/w0/0_0" in h5f
        with pytest.raises(X0SourceError):
            x0_gk.detect_representation(h5f)  # incomplete until finish()


def test_detect_saved_and_none(tmp_path):
    path = str(tmp_path / "f.h5")
    with h5py.File(path, "w") as h5f:
        assert x0_gk.detect_representation(h5f) is None
        h5f.create_group("bse/input/X0_q")
        assert x0_gk.detect_representation(h5f) == "saved"


def test_remove_x0_source(gk_file):
    with h5py.File(gk_file, "a") as h5f:
        x0_gk.remove_x0_source(h5f)
        assert "bse/x0_source" not in h5f
        assert "bse/input/X0_loc" in h5f
        x0_gk.remove_x0_source(h5f)  # no-op when absent


def test_kq_cache_is_bounded(gk_file):
    src = GkX0Source(gk_file)
    for q in np.ndindex(*DIV):
        src.get(0, "%02d.%02d.%02d" % q)
    assert len(src._kq_cache) <= x0_gk._KQ_CACHE_SIZE


def test_dcore_chiq_passes_x0q_mode_to_the_writer():
    """dcore_chiq and its worker need DCore to run; check the option plumbing statically."""
    import ast

    pkg = os.path.dirname(x0_gk.__file__)
    worker = ast.parse(open(os.path.join(pkg, "dcore_chiq_worker.py")).read())
    calls = [n for n in ast.walk(worker) if isinstance(n, ast.Call)
             and getattr(n.func, "attr", None) == "save_X0q_for_bse"]
    assert len(calls) == 1
    kw = {k.arg: ast.dump(k.value) for k in calls[0].keywords}
    assert "x0q_mode" in kw and "'X0q_mode'" in kw["x0q_mode"]

    cli = open(os.path.join(pkg, "cli", "dcore_chiq.py")).read()
    assert "params['X0q_mode'] = self._params['bse']['X0q_mode']" in cli
    assert 'pars.add_option("bse", "X0q_mode", str, "save",' in cli
