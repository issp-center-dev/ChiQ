# Design: optional on-the-fly X0(q) reconstruction from G(k)

**Date:** 2026-10-11
**Status:** Ready for human review (design review rounds 1-3 addressed)
**Base branch:** `design/pip-packaging` (PR #7)

## 1. Problem

The intermediate file `dmft_bse.h5` written by `dcore_chiq` is dominated by
`bse/input/X0_q`, which is stored densely for every saved q and every bosonic
frequency. As the number of orbitals grows, this becomes the storage bottleneck.

### 1.1 Notation (used throughout)

| Symbol | Meaning |
|---|---|
| `N_in` | inner dimension of one single-particle block (`len(gf_struct[block])`); `N` orbitals without SOC, `2N` with SOC |
| `n_c` | number of correlated shells |
| `n_b1` | number of single-particle blocks (2 without SOC: up/down; 1 with SOC) |
| `num_wf` | `[bse] num_wf` = positive fermionic cutoff; X0 stores `n_wf2 = 2*num_wf` frequencies |
| `num_wb` | `[bse] num_wb`; bosonic indices `wb = 0 .. num_wb-1` |
| `n_iwG` | frequencies kept in G(k): `2*(num_wf + num_wb - 1)` |
| `div` | k mesh `(nk0, nk1, nk2)`; `n_k = nk0*nk1*nk2` |
| `n_q` | number of saved q points |
| `n_bp` | number of stored two-particle block pairs (`b1==b3 and b2==b4`); `4*n_c^2` for local order parameters without SOC, `n_c^2` with SOC |

All sizes below are **uncompressed**, complex128 (16 B).

### 1.2 Current sizes

- `X0_q`: `16 * n_bp * N_in^4 * n_wf2 * n_q * num_wb` bytes.
  Example (no SOC, `n_c=1`, `N=10`, `num_wf=20`, `num_wb=1`): 25.6 MB per q;
  7.4 GB for `quadrant` (289 q), 26 GB for all 1024 q of a 32x32 mesh.
- Proposed G(k) store: `16 * n_b1 * n_c^2 * N_in^2 * n_iwG * n_k` bytes.
  Same example: 131 MB (`num_wf=20`), 262 MB (`num_wf=40`).
  Grows with `n_c^2`; the single-shell figure is not a general guarantee.

## 2. Goal and non-goals

**Goal.** Add an opt-in mode in which `dcore_chiq` stores G(k) instead of X0(q),
and the downstream programs reconstruct X0(q) for exactly the (wb, q) they need.

**Non-goals.**
- The default (`X0q_mode = "save"`) keeps the current behavior and file layout
  byte-for-byte in what it writes. The legacy writer in `sumk_dft_chi.py` is **not
  refactored** in this change (avoids regressing the default path; see §8).
- No change to the BSE algorithm, to `X_loc`, or to the solver backends.
- `calc_Iq_scl` does not support the new mode in this change (§6.3).
- FFT-batched reconstruction for many-q runs, IR basis, complex64 storage: future work.
- Unrelated issues found during the survey (all-q buffer at `sumk_dft_chi.py:692`,
  `G_k_iw_temp.h5` never deleted in save mode) are separate follow-ups.

## 3. User-facing interface

### 3.1 `dcore_chiq` (`[bse]` section)

```
X0q_mode = "save"   # default: current behavior (store X0_q)
X0q_mode = "gk"     # store G(k) + metadata; do not store X0_q
```

- In `"gk"` mode, `X0q_qpoints_saved` is ignored (with a printed note); any q on
  the k mesh can be requested later.
- `h5_compression*` options apply to the G(k) datasets.
- `use_temp_file` keeps working: G(k) is streamed from the temp file into the
  output file, and the temp file is deleted afterwards in `"gk"` mode.

### 3.2 `chiq_main`

No new option. The input representation is detected from the file (§4.3).
In `"gk"` files, any q on the mesh listed in `omega_q` is computed.
`omega_q` is required in gk mode (error if missing). Today `omega_q = None`
already fails at `os.path.abspath(None)` (`chiq_main.py:434`), so the
"compute everything" branch is unreachable. Saved mode keeps that behavior
unchanged.
Supported calculation types: `chi0`, `bse`, `scl` (read X0), and `rpa`/`rrpa`
(read `output/chi0_q`, unchanged: run `chi0` or `bse` first as today).

### 3.3 Compatibility contract

- New readers read old (`"save"`) files unchanged.
- Old readers fail on `"gk"` files at `check_input` with
  "Subgroup 'X0_q' not found". This is documented; no transparent fallback.
- `chiq_post`, `plot_chiq_path`, `chiq_fft`, `calc_Iq` read only output datasets
  and are unaffected.

## 4. HDF5 schema

### 4.1 Location

A new sibling group `bse/<groupname>/x0_source`, next to `info`, `input` and `output`.
`h5BSE.get_keylist_data` only scans `input`/`output`, so the new group is never
mistaken for a registered quantity (`h5bse.py:346`). The registry `_keyList` is
not extended.

### 4.2 Contents

```
bse/<g>/x0_source/                   (group)
  attrs:
    format        = "gk"
    version       = 1
    complete      = 0 / 1            # set to 1 only after everything is written
    num_wf        = int              # X0 window: fermionic indices [-num_wf, num_wf)
    list_wb       = int[]            # bosonic indices supported
    n_iwG         = int              # length of G frequency axis
    iw_offset     = int              # array index of fermionic index 0 (= n_iwG // 2)
    div           = int[3]
    beta          = float
  block1_names    : str[n_b1]        # single-particle block names, in order
  block2_pairs    : int[n_b2, 4]     # (icrsh1, b1, icrsh2, b2) per two-particle block,
                                     #   same order as info/block_name
  inner2_pairs    : int[n_in2, 2]    # (i1, i2) per two-particle inner index,
                                     #   same order as info/inner_name
  G/<b1 index>    : complex128[n_c, n_c, N_in, N_in, n_iwG, n_k]
                    # layout of Gk_corr[bname] (sumk_dft_chi.py:879);
                    # chunks = (1, 1, 1, 1, n_iwG, n_k)
```

k index flattening: `k = kz + div[2]*(ky + div[1]*kx)` (`sumk_dft_chi.py:341`).

### 4.3 Representation detection and stale data

The reader decides with this table (all checks before any rank-dependent work):

| `input/X0_q` | `x0_source` | result |
|---|---|---|
| present | absent | saved mode (current code path) |
| absent | present, `complete=1`, `format="gk"`, `version=1` | gk mode |
| absent | present but incomplete / unknown format or version | error |
| present | present | error: ambiguous representation; tells the user to rerun `dcore_chiq` |
| absent | absent | current error ("Subgroup 'X0_q' not found") |

Writer policy (`dcore_chiq`, rank 0):
- `"gk"` run: delete `input/X0_q` and any old `x0_source`, then write
  `x0_source` with `complete=0`, write all data, and set `complete=1` last.
  An interrupted run therefore leaves an incomplete group and is rejected.
- `"save"` run: delete any existing `x0_source` before writing X0_q.
- Both: print a note that HDF5 does not reclaim freed space and suggest
  `h5repack` or a fresh output file to actually shrink an existing file.

Other `input/*` quantities (`X_loc`, `chi_loc_in`, `gamma0`, ...) are untouched.

Validation by `GkX0Source` when opening (any failure is an error naming the check):
- `G/<b>` exists for every `b < len(block1_names)`, with shape
  `(n_c, n_c, N_in, N_in, n_iwG, prod(div))`, complex dtype.
- `block2_pairs` is an integer array of shape `(n_b2, 4)` with
  `0 <= icrsh < n_c` and `0 <= b < len(block1_names)`; `n_b2 == len(info/block_name)`.
- `inner2_pairs` is an integer array of shape `(n_in2, 2)` with
  `0 <= i < N_in`; `n_in2 == len(info/inner_name)`. (Negative values are rejected
  explicitly: NumPy would silently accept them as indices from the end.)
- `beta` equals `info/beta`.
- `iw_offset == n_iwG // 2`, and for every `wb` in `list_wb`:
  `0 <= wb` and `iw_offset + num_wf + wb <= n_iwG` and `iw_offset - num_wf >= 0`.
- `X0_loc` exists for every `wb` in `list_wb` (needed by `_get_niw`).

**Outputs are not versioned against inputs.** `rpa`/`rrpa` read `output/chi0_q`
from the output file. After rerunning `dcore_chiq` (in either mode), the user
must recompute `chi0`/`bse` into a fresh output file. Otherwise old and new
results can mix, as is already the case today. This is documented in the
`X0q_mode` option description and in the docs; no provenance check is added.

## 5. Mathematics (the contract the reconstruction must satisfy)

For a stored two-particle block pair `(b_l, b_r)` with
`block2_pairs[b_l] = (c1, b1, c2, b2)`, `block2_pairs[b_r] = (c3, b3, c4, b4)`,
kept only if `b1 == b3 and b2 == b4` (`sumk_dft_chi.py:661,700`), inner indices
`i = (i1, i2)`, `j = (i3, i4)`, bosonic index `wb >= 0`, and X0 frequency
position `t = 0 .. 2*num_wf-1`:

```
gl(t, k) = G[b1][c3, c1, i3, i1, iw_offset - num_wf + t,      k]
gr(t, k) = G[b2][c2, c4, i2, i4, iw_offset - num_wf + t + wb, k]

X0_q[(b_l, b_r)][i, j, t] = -(1/n_k) * sum_k gl(t, k) * gr(t, (k-q) mod div)
```

- Minus sign, no complex conjugation, `1/n_k` normalization, no `1/beta`;
  output axis order `(i, j, t)`.
- **q sign.** The formula uses `k-q` to match, up to rounding, the default FFT
  writer (`calc_chi0_fft`, `sumk_dft_chi.py:403`). That writer computes
  `fft( ifft(gl)(r) * ifft(gr)(-r) )`, which equals `sum_k gl(k) gr(k-q)`; this was
  verified numerically. The legacy direct-sum writer (`calc_chi0_sum`) uses `k+q`.
  It is not reachable from `dcore_chiq` (the worker always uses `algo="fft"`).
  The two agree only when X0(q) = X0(-q), e.g. with inversion symmetry. The
  discrepancy in the unused legacy path is reported separately and not changed here.
- `wb` must be in `list_wb`; then the slice stays inside `[0, n_iwG)`.
- The output has the full `2*num_wf` window. `chiq_main` keeps truncating with
  its current leading slice `[..., :n_iw]` (`chiq_main.py:215,227`), so
  `[chiq_main] num_wf` behaves exactly as with saved files.

**q labels.** `omega_q` labels are `"qx.qy.qz"`. In gk mode they are parsed as
integers with `0 <= q_l < div[l]`. Labels off the mesh are rejected before any
work. Output datasets keep the label string as given, so `chiq_post` keys match.

**X0_loc in gk mode.** The current FFT writer stores the average of X0 over **all**
q of the mesh (`sumk_dft_chi.py:716`). Over the full mesh this factorizes exactly:

```
X0_loc[(b_l, b_r)][i, j, t] = -( mean_k gl(t, k) ) * ( mean_k gr(t, k) )
```

`dcore_chiq` computes it this way and stores it in `input/X0_loc` as today.
The k-averages `Gbar[b][c, c', i, i', w] = mean_k G[b][c, c', i, i', w, k]` are
accumulated while G is written to `x0_source` (§6.2). This works both for the
in-memory `Gk_corr` and for the `use_temp_file` path, where `Gk_corr` is `None`
and G exists only in the temp file. `chi0_loc_in` keeps its current computation from the
full-window `G_loc_corr` (`sumk_dft_chi.py:791`); it must not be derived from the
short window.

## 6. Code structure

### 6.1 New DCore-free module `chiq/x0_gk.py`

- `parse_q_label(label, div) -> (qx, qy, qz)` with validation.
- `kq_permutation(div, q) -> int[n_k]`, vectorized (`np.indices` + modular
  arithmetic), computed once per q.
- `x0_block(Gb1, Gb2, pair_l, pair_r, inner2_pairs, perm, wb, num_wf, iw_offset)
  -> complex[n_in2, n_in2, 2*num_wf]`, implemented as one batched contraction
  per frequency: `(n_in2 x n_k) @ (n_k x n_in2)` GEMMs, without Python loops over
  orbitals or k.
- `x0_loc_block(...)` for the factorized local average.
- `class GkX0Source`: opens and validates `x0_source` (§4.3), loads G into memory
  once per process, and exposes
  - `wb_list()`,
  - `get(wb, q_label) -> {(b_l, b_r): ndarray}`, in the same h5bse format as
    `h5BSE.get(("X0_q", wb, q))`, so `_convert_block_matrix` is reused unchanged.

### 6.2 `dcore_chiq` / `sumk_dft_chi.py`

- `dcore_chiq.py`: new option, passed as `params['X0q_mode']` to the worker.
- `SumkDFTChi.save_X0q_for_bse(..., x0q_mode="save")`: in `"gk"` mode, after
  `extract_Gk_crsh`:
  1. `_save_x0_source(temp_file_ar)` (rank 0) writes §4.2. For each
     `(bname, c1, c2, in1, in2)` it takes the `(n_iwG, n_k)` slice from
     `Gk_corr` or from the temp file (`ds_name` of `sumk_dft_chi.py:968`), writes
     it into the `G/<b>` dataset, and stores its k-mean into `Gbar`.
     Peak extra memory is one slice plus `Gbar`.
  2. Per wb: `_save_X0_loc_from_gbar(wb)` using `Gbar`, then `_save_chi0_loc(wb)`
     as today.
  3. Set `complete=1`, close the archive, and only then delete the temp file.
  `_save_X0_data` is not called. The `"save"` branch is not changed.

### 6.3 `chiq_main` and other consumers

- `chi_base` gets an X0 provider chosen once at start (§4.3): `SavedX0Source`
  (wraps the existing `h5BSE.get`) or `GkX0Source`.
- `check_input` checks provider availability instead of `X0_q in datalist` for
  types whose `input_q == "X0_q"`.
- **Target scheduling.** The main loop iterates exactly over this rank's
  `(omega, q)` tuples. Today it filters omega and q independently, so a rank
  can compute combinations it did not get (`chiq_main.py:537-566`). In saved mode,
  tuples missing from the file are skipped with a log line, as today. This fix
  applies to both modes; it is a separate commit with its own test.
- The reconstructed X0 is a local variable per (omega, q); it is not retained.
- `calc_Iq_scl`: if `x0_source` is present, exit with an error saying that gk
  files are not supported by this command yet and to use `X0q_mode = "save"`.

### 6.4 Memory per process (gk mode)

`G store (§1.2) + n_bp * 16 * N_in^4 * 2*num_wf (one X0(q)) + solver matrices
(unchanged) + accumulated outputs (unchanged)`. G is replicated on every rank;
this is documented together with the formula. Lazy per-block loading is future
work if this becomes a limit.

## 7. Performance

- Direct sum per q: `n_bp * n_in2^2 * 2*num_wf * n_k` complex multiply-adds,
  executed as GEMMs. FFT over all q would be `O(n_k log n_k)` per orbital
  combination. Direct sum wins when `n_q << log2(n_k) * const`; for sparse q
  sets (paths) it is the right choice.
- Acceptance: a benchmark script (not in CI) on synthetic data
  (`N_in=10`, `n_k=1024`, `num_wf=20`) reports, per q, reconstruction time vs.
  reading the same X0 from a saved file, and vs. the `bse` solve time.
  The numbers go in the PR description. "Negligible" is not claimed without them.

## 8. Testing

Tolerances: raw X0 and X0_loc `rtol=1e-12`; final observables `rtol=1e-10`
(FFT vs. direct-sum rounding can be amplified near poles).

**DCore-free unit tests** (`tests/python/non-mpi/x0_gk/`):
1. Kernel vs. an independent oracle: explicit Python loops over k using `k-q`,
   written separately from the kernel. Includes a case where X0(q) != X0(-q)
   to pin the sign. Random complex G with asymmetric off-diagonal
   elements; non-cubic mesh `div=(4,3,2)`; `n_c=2`; a SOC-like single block with
   `N_in=4`; `wb = 0` and `wb = max`; nonzero q; nonlocal block pairs (`c1 != c2`).
2. Kernel vs. an independent FFT formula (ifft/fft with the `-r` inversion of
   `sumk_dft_chi.py:403`) on the same G. This pins the convention of the legacy
   FFT writer.
3. `X0_loc` factorized formula vs. the mean over all q of the oracle.
4. `GkX0Source` round trip through a synthetic file, plus rejection cases:
   incomplete group, unknown version, both representations present, q off the
   mesh, wb not in `list_wb`, malformed mappings (wrong shape or dtype, negative or
   out-of-range shell/block/inner indices).
5. `chiq_main` end-to-end on synthetic files: a gk file and an equivalent saved
   file (X0_q materialized from the oracle) give identical outputs for `chi0`,
   `bse`, `scl`, then `rpa` and `rrpa`. Includes `[chiq_main] num_wf` truncation.
6. Scheduling: a rank does not compute tuples it was not assigned.
7. `calc_Iq_scl` on a gk file exits with the documented message.

**MPI** (`tests/python/mpi/`): test 5 with 2 ranks, gk mode.

**DCore integration** (CI scientific job, where DCore is installed): run the
X0 part of `dcore_chiq` (`skip_Xloc = true`) in both modes on a small example.
Then compare saved `X0_q` with the `GkX0Source` reconstruction for all saved q,
and `X0_loc` / `chi0_loc_in` between modes. Run it with `use_temp_file = false`
and `true`, and check that the temp file is removed in gk mode. This is the
regression guard for the untouched legacy writer.

## 9. Review log

- Round 1 (design review): schema and discovery, frequency window, block/q/SOC
  mapping, X0_loc/chi0_loc_in provenance, consumer coverage, stale-data policy,
  test oracle, cost claim, MPI scheduling, size definitions. All addressed above.
- Round 2: q sign (FFT writer is `k-q`, verified numerically; §5 now uses `k-q`),
  temp-file path for X0_loc (§5, §6.2), explicit validation list (§4.3),
  `omega_q` required in gk mode (§3.2), output staleness documented (§4.3).
- Round 3: no blocking issues; implementation-ready. Added dtype, shape and
  non-negativity checks for the index mappings and matching rejection tests (§4.3, §8).
