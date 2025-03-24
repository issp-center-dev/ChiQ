# Original script: Ce_gamma.py from
# https://triqs.ipht.cnrs.fr/1.4/applications/dft_tools/guide/images_scripts/Ce-gamma_script.html

import sys
import os

import numpy as np

from chiq.sumk_dft_chi import *
from dcorelib.triqs_compat.gf import BlockGf, GfImFreq

#=============================================================================
# set input file

argvs = sys.argv
argc = len(argvs)
if argc != 2:
    raise Exception("How to use: python3 %s file_param" %argvs[0])

file_param=argvs[1]
if not os.path.exists(file_param):
    raise Exception("File '%s' not exists" %file_param)

#=============================================================================
# set paramters

def str2bool(str):
    return str in ['True', 'true', '1', 'y', 'Y']

import configparser

default_config = {
    # 'dft_h5file' : 'hubbard.h5',
    'n_wf_G1' : '1024'
}
config = configparser.ConfigParser(default_config)
config.read(file_param)

# DMFT params
beta = float(config.get('DMFT', 'beta'))
U_int = float(config.get('DMFT', 'U_int'))
dft_h5file = config.get('DMFT', 'dft_h5file')
n_wf_G1 = int(config.get('DMFT', 'n_wf_G1'))

# After DMFT loop
n_wf_G2 = int(config.get('POSTPROCESS', 'n_wf_G2'))

# fixed parameters
J_hund = 0.0
chemical_potential_init=0
fix_mu = True

del config

if mpi.is_master_node():
    print("file_param =", file_param)
    print("[DMFT]")
    print(" beta =", beta)
    print(" U_int =", U_int)
    print(" dft_h5file =", dft_h5file)
    print(" n_wf_G1 =", n_wf_G1)
    print("[POSTPROCESS]")
    print(" n_wf_G2 =", n_wf_G2)

#=============================================================================

# Init the SumK class
SK=SumkDFT(hdf_file=dft_h5file, use_dft_blocks=False)
# dft_data==dft_data_fbz since IBZ is not used

Norb = SK.corr_shells[0]['dim']
l    = SK.corr_shells[0]['l']

# gf_struct = { sn: range(-l,l+1) for sn in ['up', 'down'] }
gf_struct = SK.gf_struct_solver[0]
spin_names = list(gf_struct.keys())
orb_names = gf_struct[spin_names[0]]

# Init Sigma_iw
block_list = [GfImFreq(indices=indices, beta=beta, n_points=n_wf_G1) for name, indices in gf_struct.items()]
Sigma_iw = BlockGf(name_list=spin_names, block_list=block_list, make_copies=True)
Sigma_iw.zero()

# Hamiltonian
# U_mat = U_matrix(l=l, U_int=U_int, J_hund=J_hund, basis='spherical')
U_mat = np.array([[[[U_int]]]])
# h_int = hamiltonian.h_int_slater(spin_names, orb_names, U_matrix=U_mat, off_diag=True)

chemical_potential=chemical_potential_init
# load previous data: old self-energy, chemical potential, DC correction

if fix_mu:
    SK.set_mu(chemical_potential)
else:
    chemical_potential = SK.calc_mu( precision = 0.000001 )

SK.set_Sigma([ Sigma_iw ])  # for Loops=0

#--------------------------------------------------
# chi0_q

if True:
    SKC = SumkDFTChi(hdf_file=dft_h5file, dft_data_fbz='dft_input', use_dft_blocks=False)
    SKC.set_mu(SK.chemical_potential)
    SKC.set_dc(SK.dc_imp, SK.dc_energ)
    SKC.set_Sigma([ Sigma_iw ])

    # SK.init_for_X0()
    wbs=[0, 3]

    # FFT
    SKC.save_X0q_for_bse(list_wb=wbs, n_wf_cutoff=n_wf_G2, h5_file='chi0q_fft.h5', flag_save_X0_loc=False, flag_save_chi0_loc=False, qpoints_saved='quadrant', del_hopping=False)

    # sum
    q_dict = {}
    for x, y in product(list(range(3)), list(range(3))):
        q_label = "%02d.%02d.00" % (x, y)
        q_dict[q_label] = (x,y,0)

    SKC.save_X0q_for_bse(list_wb=wbs, n_wf_cutoff=n_wf_G2, h5_file='chi0q_sum.h5', algo='sum', q_dict=q_dict, flag_save_X0_loc=False, flag_save_chi0_loc=False)

#--------------------------------------------------
