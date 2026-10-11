import os

import h5py
import numpy as np

from chiq import h5bse
from chiq.x0_gk import GkX0Source, detect_representation

saved = h5bse.h5BSE('x0_save.h5')
keys = saved.get_keylist_data('input')
x0q = [k for k in keys if k[0] == 'X0_q']
assert x0q, "no X0_q in x0_save.h5"

for gk_file in ('x0_gk.h5', 'x0_gk_temp.h5'):
    with h5py.File(gk_file, 'r') as f:
        assert detect_representation(f) == 'gk'
        assert 'bse/input/X0_q' not in f
    src = GkX0Source(gk_file)
    gk = h5bse.h5BSE(gk_file)
    for _, wb, q in x0q:
        ref = saved.get(('X0_q', wb, q))
        got = src.get(wb, q)
        assert set(ref) == set(got), (wb, q)
        for blk in ref:
            np.testing.assert_allclose(got[blk], ref[blk], rtol=1e-12, atol=1e-14, err_msg=f"{gk_file} X0_q {wb} {q} {blk}")
    for name in ('X0_loc', 'chi0_loc_in'):
        for wb in (0, 3):
            ref = saved.get((name, wb))
            got = gk.get((name, wb))
            assert set(ref) == set(got)
            for blk in ref:
                np.testing.assert_allclose(got[blk], ref[blk], rtol=1e-12, atol=1e-14, err_msg=f"{gk_file} {name} {wb} {blk}")

assert not os.path.exists('G_k_iw_temp.h5'), "temp file must be removed in gk mode"
print("X0 MODES AGREE")
