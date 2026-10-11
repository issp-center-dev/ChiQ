"""Write X0 with X0q_mode='save' and 'gk' (with and without temp file) from the same model."""

import os
import sys

from chiq.sumk_dft_chi import SumkDFTChi, SumkDFT, mpi
from dcorelib.triqs_compat.gf import BlockGf, GfImFreq

import configparser

config = configparser.ConfigParser({'n_wf_G1': '1024'})
config.read(sys.argv[1])
beta = float(config.get('DMFT', 'beta'))
dft_h5file = config.get('DMFT', 'dft_h5file')
n_wf_G1 = int(config.get('DMFT', 'n_wf_G1'))
n_wf_G2 = int(config.get('POSTPROCESS', 'n_wf_G2'))

SK = SumkDFT(hdf_file=dft_h5file, use_dft_blocks=False)
gf_struct = SK.gf_struct_solver[0]
spin_names = list(gf_struct.keys())
block_list = [GfImFreq(indices=indices, beta=beta, n_points=n_wf_G1) for name, indices in gf_struct.items()]
Sigma_iw = BlockGf(name_list=spin_names, block_list=block_list, make_copies=True)
Sigma_iw.zero()
SK.set_mu(0)
SK.set_Sigma([Sigma_iw])

wbs = [0, 3]
runs = [
    ('x0_save.h5', 'save', None),
    ('x0_gk.h5', 'gk', None),
    ('x0_gk_temp.h5', 'gk', 'G_k_iw_temp.h5'),
]
for h5_file, mode, temp_file in runs:
    # a fresh instance per run: save_X0q_for_bse reduces Sigma in place
    SKC = SumkDFTChi(hdf_file=dft_h5file, dft_data_fbz='dft_input', use_dft_blocks=False)
    SKC.set_mu(SK.chemical_potential)
    SKC.set_dc(SK.dc_imp, SK.dc_energ)
    SKC.set_Sigma([Sigma_iw])
    SKC.save_X0q_for_bse(list_wb=wbs, n_wf_cutoff=n_wf_G2, h5_file=h5_file, qpoints_saved='quadrant',
                         temp_file=temp_file, x0q_mode=mode)
