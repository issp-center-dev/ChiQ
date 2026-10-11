"""chiq_main with 2 MPI ranks on a G(k) file computes exactly the requested targets."""

import os
import subprocess
import sys

import h5py
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "non-mpi", "x0_gk"))
import synth_gk as S  # noqa: E402

import chiq  # noqa: E402
from chiq import h5bse  # noqa: E402

PACKAGE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(chiq.__file__)))
# Sorted and split over 2 ranks as [3, 2] (more_itertools.divide), so rank 0 gets both
# frequencies with different q: independent omega/q filtering would add (1, 00.00.00),
# (1, 01.00.00) and (0, 02.01.00) on rank 0.
TARGETS = [("0", "00.00.00"), ("0", "01.00.00"), ("1", "02.01.00"), ("1", "00.01.00"), ("1", "01.01.00")]


def _toml(path, infile, outfile):
    path.write_text('[chiq_common]\ninput = "%s"\noutput = "%s"\ntype = ["chi0"]\nomega_q = "q.dat"\n'
                    '[chiq_main]\nbackend = "numpy"\n' % (infile, outfile))
    return str(path)


def _run(nproc, toml, cwd):
    env = dict(os.environ)
    env["PYTHONPATH"] = PACKAGE_DIR + os.pathsep + env.get("PYTHONPATH", "")
    env.setdefault("OMPI_MCA_rmaps_base_oversubscribe", "1")
    cmd = ["mpirun", "-np", str(nproc), sys.executable, "-m", "chiq.cli.chiq_main", toml]
    return subprocess.run(cmd, cwd=cwd, env=env, capture_output=True, text=True)


def _computed_per_rank(cwd, nproc):
    """(omega, q) computed by each rank, from the chiq_main_<rank>.log files."""
    per_rank = []
    for rank in range(nproc):
        done, omega = [], None
        with open(os.path.join(cwd, "chiq_main_%d.log" % rank)) as f:
            for line in f:
                msg = line.rsplit(":", 1)[-1].strip()
                if msg.startswith("omega = "):
                    omega = msg.split("=")[1].strip()
                elif msg.startswith("q = "):
                    done.append((omega, msg.split("=")[1].strip()))
        per_rank.append(done)
    return per_rank


def _chi0q(path):
    with h5py.File(path, "r") as f:
        g = f["bse/output/chi0_q"]
        return {(w, q): {b: g[w][q][b][()] for b in g[w][q]} for w in g for q in g[w]}


def test_two_ranks_compute_exactly_the_targets(tmp_path):
    rng = np.random.default_rng(11)
    div, num_wf, list_wb = (3, 2, 1), 2, [0, 1]
    n_iwG = 2 * (num_wf + 1)
    G = S.random_g(rng, 2, 1, 1, n_iwG, 6)
    block2, inner2 = S.make_pairs(1, 2, 1, nonlocal_pairs=False)
    S.write_gk_file(str(tmp_path / "gk.h5"), G, block2, inner2, num_wf, list_wb, div, 5.0)

    # equivalent saved file (X0_q from the independent FFT formula)
    S.write_info(str(tmp_path / "saved.h5"), [str(p) for p in block2], [str(p) for p in inner2], 5.0)
    bs = h5bse.h5BSE(str(tmp_path / "saved.h5"))
    for wb in list_wb:
        fft = S.fft_x0_all_q(G, block2, inner2, wb, num_wf, n_iwG // 2, div)
        for iq, q in enumerate(np.ndindex(*div)):
            bs.save(key=("X0_q", wb, "%02d.%02d.%02d" % q), data={k: v[..., iq] for k, v in fft.items()})
        bs.save(key=("X0_loc", wb), data={k: v.mean(axis=-1) for k, v in fft.items()})

    (tmp_path / "q.dat").write_text("".join("%s %s\n" % t for t in TARGETS))
    r = _run(2, _toml(tmp_path / "gk.toml", "gk.h5", "gk.out.h5"), str(tmp_path))
    assert r.returncode == 0, r.stdout + r.stderr
    per_rank = _computed_per_rank(str(tmp_path), 2)  # before the next run overwrites the logs
    r = _run(1, _toml(tmp_path / "saved.toml", "saved.h5", "saved.out.h5"), str(tmp_path))
    assert r.returncode == 0, r.stdout + r.stderr

    assert sorted(per_rank[0] + per_rank[1]) == sorted(TARGETS)  # no extras, no duplicates
    assert len(per_rank[0]) == 3 and len(per_rank[1]) == 2
    assert {w for w, _ in per_rank[0]} == {"0", "1"}  # the mixed-frequency rank

    gk = _chi0q(str(tmp_path / "gk.out.h5"))
    saved = _chi0q(str(tmp_path / "saved.out.h5"))
    assert set(gk) == {("w" + w, "q_" + q) for w, q in TARGETS}
    for key, blocks in gk.items():
        for b, value in blocks.items():
            np.testing.assert_allclose(value, saved[key][b], rtol=1e-10, atol=1e-13)
