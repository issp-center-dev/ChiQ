# from triqs_dft_tools.sumk_dft import SumkDFT
import sys
import os
import time
import pprint
import numpy as np
from itertools import product
from mpi4py import MPI
from .index_pair import IndexPair, IndexPair2
from .h5bse import h5BSE
from .tools import WallTime

# from .pytriqs_gf_compat import *
# from .dft_tools_compat import SumkDFT
from dcore._dispatcher import HDFArchive, mpi, BlockGf, GfImFreq, MeshImFreq
from dcore.sumkdft_opt import SumkDFT_opt as SumkDFT

from dcore.mpi import gatherv


def get_freq_index(n_iw, iw_cutoff=(None, None)):
    """
    return:
        (int)wmin, (int)wmax: Indices (positions in array) of the range of frequencies

        Ex.) If iw_cutoff is not given, (w_min, wmax) = (0, n_iw)

    input:
        (int)n_iw:  Number of frequencies including negative side
        (int, int)iw_cutoff:  Lower and upeer frequency cutoff (Matsubara index)
    """

    assert n_iw % 2 == 0  # even number
    w0 = n_iw // 2
    wmin = 0 if iw_cutoff[0] == None else w0 + iw_cutoff[0]
    wmax = 2 * w0 if iw_cutoff[1] == None else w0 + iw_cutoff[1]
    assert wmin >= 0
    assert wmax <= 2 * w0
    return wmin, wmax


def get_block_size(g):
    """
    Get block size of g

    Parameters
    ----------
    g : object of GfImFreq, GfImTime or GfLegendre

    Returns
    -------
    block_size : int
        block of g (e.g. number of orbitals)
    """

    assert len(g.indices[0]) == len(g.indices[1])
    return len(g.indices[0])


def gf_reduce_iw(gf_iw, n_points):
    """
    Reduce number of Matsubara points

    gf_iw  : BlockGf
    return : BlockGf
    """

    assert isinstance(gf_iw, BlockGf), f"type(gf_iw)={type(gf_iw)}"
    for bname, g in gf_iw:
        assert isinstance(g, GfImFreq), f"type(g)={type(g)}"

    # creat a new mesh
    reduced_mesh = GfImFreq(indices=[0, ], beta=gf_iw.mesh.beta, n_points=n_points).mesh

    def is_equal_mesh_except_npoints(mesh1, mesh2):
        """check if two mesh are equivalent excepting 'n_points'"""
        return all([mesh1.beta == mesh2.beta, mesh1.statistic == mesh2.statistic])
        # positive_only is introduced in ver1.5

    assert is_equal_mesh_except_npoints(reduced_mesh, gf_iw.mesh)

    gf_struct = {bname: g.indices for bname, g in gf_iw}
    gf_iw_reduced = BlockGf(name_block_generator=[(block, GfImFreq(indices=inner, mesh=reduced_mesh))
                                                  for block, inner in list(gf_struct.items()) ],
                            make_copies=False)

    # range of iw to be extracted
    bname0 = list(gf_iw.indices)[0]  # the first block name
    wmin, wmax = get_freq_index(n_iw=gf_iw[bname0].data.shape[0], iw_cutoff=(-n_points, n_points))

    # copy data in the reduced range
    for bname, g in gf_iw:
        gf_iw_reduced[bname].data[:, :, :] = g.data[wmin:wmax, :, :]
    return gf_iw_reduced


class SumkDFTChi_aux(SumkDFT):
    """
    Auxiliary class for SumkDFTChi
    """

    # __init__ is common with the super class

    def downfold_offdiagonal(self, ik, ish1, ish2, bname, gf_to_downfold, gf_inp, shells='corr', ir=None):
        r"""
        [ORIGINAL: function 'downfold' in sumk_dft.py]
        Downfolds a block of the Green's function for a given shell and k-point using the corresponding projector matrices.

        Parameters
        ----------
        ik : integer
             k-point index for which the downfolding is to be done.
        ish : integer
              Shell index of GF to be downfolded.

              - if shells='corr': ish labels all correlated shells (equivalent or not)
              - if shells='all': ish labels only representative (inequivalent) non-correlated shells

        bname : string
                Block name of the target block of the lattice Green's function.
        gf_to_downfold : Gf
                       Block of the Green's function that is to be downfolded.
        gf_inp : Gf
                 FIXME
        shells : string, optional

                 - if shells='corr': orthonormalized projectors for correlated shells are used for the downfolding.
                 - if shells='all': non-normalized projectors for all included shells are used for the downfolding.

        ir : integer, optional
             Index of equivalent site in the non-correlated shell 'ish', only used if shells='all'.

        Returns
        -------
        gf_downfolded : Gf
                      Downfolded block of the lattice Green's function.
        """

        gf_downfolded = gf_inp.copy()
        # get spin index for proj. matrices
        isp = self.spin_names_to_ind[self.SO][bname]
        n_orb = self.n_orbitals[ik, isp]
        if shells == 'corr':
            dim1 = self.corr_shells[ish1]['dim']
            dim2 = self.corr_shells[ish2]['dim']
            projmat1 = self.proj_mat[ik, isp, ish1, 0:dim1, 0:n_orb]
            projmat2 = self.proj_mat[ik, isp, ish2, 0:dim2, 0:n_orb]
        elif shells == 'all':
            if ir is None:
                raise ValueError("downfold: provide ir if treating all shells.")
            dim1 = self.shells[ish1]['dim']
            dim2 = self.shells[ish2]['dim']
            projmat1 = self.proj_mat_all[ik, isp, ish1, ir, 0:dim1, 0:n_orb]
            projmat2 = self.proj_mat_all[ik, isp, ish2, ir, 0:dim2, 0:n_orb]

        gf_downfolded.from_L_G_R(
            projmat1, gf_to_downfold, projmat2.conjugate().transpose())

        return gf_downfolded

    def extract_G_loc_offdiagonal(self, mu=None, iw_or_w='iw', with_Sigma=True, with_dc=True, broadening=None):
        r"""
        [ORIGINAL: function 'extract_G_loc' in sumk_dft.py]
        Extracts the local downfolded Green function by the Brillouin-zone integration of the lattice Green's function.

        Parameters
        ----------
        mu : real, optional
             Input chemical potential. If not provided the value of self.chemical_potential is used as mu.
        with_Sigma : boolean, optional
                     If True then the local GF is calculated with the self-energy self.Sigma_imp.
        with_dc : boolean, optional
                  If True then the double-counting correction is subtracted from the self-energy in calculating the GF.
        broadening : float, optional
                     Imaginary shift for the axis along which the real-axis GF is calculated.
                     If not provided, broadening will be set to double of the distance between mesh points in 'mesh'.
                     Only relevant for real-frequency GF.

        Returns
        -------
        G_loc : 2-dimensional list of BlockGf (Green's function) objects
                List of the local Green's functions G_loc[icrsh1][icrsh2]
                rotated into the corresponding local frames.

        """

        beta = self.Sigma_imp_iw[0].mesh.beta
        bname0 = list(self.Sigma_imp_iw[0].indices)[0]

        # Sigma_imp_iw is defined for correlated shells
        assert len(self.Sigma_imp_iw) == self.n_corr_shells

        # All correlated shells have the same block/inner structure
        for icrsh in range(self.n_corr_shells):
            assert self.Sigma_imp_iw[0].mesh == self.Sigma_imp_iw[icrsh].mesh

        if mu is None:
            mu = self.chemical_potential

        if iw_or_w == "iw":
            G_loc_mat = [[self.Sigma_imp_iw[0].copy() for icrsh2 in range(self.n_corr_shells)]
                     for icrsh1 in range(self.n_corr_shells)]   # this list will be returned
            # G_loc_inequiv = [BlockGf(name_block_generator=[(block, GfImFreq(indices=inner, mesh=G_loc[0].mesh)) for block, inner in self.gf_struct_solver[ish].iteritems()],
            #                          make_copies=False) for ish in range(self.n_inequiv_shells)]
        elif iw_or_w == "w":
            raise Exception("iw_or_w=='w' not implemented")

        for icrsh1, icrsh2 in product(list(range(self.n_corr_shells)), repeat=2):
            G_loc_mat[icrsh1][icrsh2].zero()                          # initialize to zero

        ikarray = np.array(list(range(self.n_k)))
        for ik in mpi.slice_array(ikarray):
            if iw_or_w == 'iw':
                G_latt = self.lattice_gf(
                    ik=ik, mu=mu, iw_or_w=iw_or_w, with_Sigma=with_Sigma, with_dc=with_dc, beta=beta)
            # elif iw_or_w == 'w':
            G_latt *= self.bz_weights[ik]

            for icrsh1, icrsh2 in product(list(range(self.n_corr_shells)), repeat=2):
                # init temporary storage
                tmp = G_loc_mat[icrsh1][icrsh2].copy()
                for bname, gf in tmp:
                    tmp[bname] << self.downfold_offdiagonal(
                        ik, icrsh1, icrsh2, bname, G_latt[bname], gf)
                G_loc_mat[icrsh1][icrsh2] += tmp

        # Collect data from mpi
        for icrsh1, icrsh2 in product(list(range(self.n_corr_shells)), repeat=2):
            G_loc_mat[icrsh1][icrsh2] << mpi.all_reduce(
                mpi.world, G_loc_mat[icrsh1][icrsh2], lambda x, y: x + y)
        mpi.barrier()

        # G_loc[:] is now the sum over k projected to the local orbitals.
        # here comes the symmetrisation, if needed:
        if self.symm_op != 0:
            for icrsh1 in range(self.n_corr_shells):
                G_loc_mat[icrsh1] = self.symmcorr.symmetrize(G_loc_mat[icrsh1])

        # G_loc is rotated to the local coordinate system:
        if self.use_rotations:
            # for r in self.rot_mat:
            #     if not np.allclose(np.identity(r.shape[0]), r):
            #         raise Exception("Not implemented: rot_mat should an identity matrix!")
            # FIXME
            for icrsh1, icrsh2 in product(list(range(self.n_corr_shells)), repeat=2):
                for bname, gf in G_loc_mat[icrsh1][icrsh2]:
                    gf.from_L_G_R(self.rot_mat[icrsh1].conjugate().transpose(), gf, self.rot_mat[icrsh2])

        return G_loc_mat


class SumkDFTChi(SumkDFTChi_aux):
    """
    Extends the SumkDFT class to susceptibility calculations.
    """

    def __init__(
        self,
        hdf_file,
        h_field=0.0,
        use_dft_blocks=False,
        dft_data="dft_input",
        symmcorr_data="dft_symmcorr_input",
        parproj_data="dft_parproj_input",
        symmpar_data="dft_symmpar_input",
        bands_data="dft_bands_input",
        transp_data="dft_transp_input",
        misc_data="dft_misc_input",
        dft_data_fbz="dft_input_fbz",
        dft_chi_data="dft_input_chi",
    ):
        """
        Initialization of the class.

        There are 2 additional parameters compared with SumkDFT class:
            dft_data_fbz : Name of hdf5 subgroup in which DFT data prepared for the full BZ (FBZ) are stored.
            chi_data : Name of hdf5 subgroup in which auxiliary data necessary of chi0 calculation are stored.

        dft_data can be either IBZ or FBZ.
        """
        args = {'hdf_file' : hdf_file,
                'h_field' : h_field,
                'use_dft_blocks' : use_dft_blocks,
                # 'dft_data' : dft_data,
                'symmcorr_data' : symmcorr_data,
                'parproj_data' : parproj_data,
                'symmpar_data' : symmpar_data,
                'bands_data' : bands_data,
                'transp_data' : transp_data,
                'misc_data' : misc_data,
                }

        # init with FBZ database
        super().__init__(dft_data=dft_data_fbz, **args)

        assert self.symm_op == 0

        # read supplemental data from subgroup dft_chi_data in HDF5
        subgrp = dft_chi_data
        with HDFArchive(hdf_file, 'r') as ar:
            self.div = ar[subgrp]['div']
            self.nk_all = np.prod(self.div)

        assert self.n_k == self.nk_all, "n_k=%d in 'dft_data' is not consistent with 'dft_chi_data' (nk_all=%d). 'dft_data' may not be in FBZ" %(self.n_k, self.nk_all)

        # make an instance of the auxiliary class using the IBZ database
        # This will be used for computing G_loc[icrsh1][icrsh2]
        if dft_data == dft_data_fbz:
            # skip making a new instance if IBZ data is not provided
            self.SK_ibz = self
        else:
            self.SK_ibz = SumkDFTChi_aux(dft_data=dft_data, **args)

        self.minus_r = self.__calc_minus_r()

        # For performance measurement
        self.__to_save_time = False
        self.__time_starts = {}
        self.__time_elapsed = {}

    def __del__(self):
        if self.__to_save_time and mpi.is_master_node():
            with open("time.dat", "w") as f:
                for key in self.__time_elapsed:
                    f.write(f"{key}: {self.__time_elapsed[key]} sec\n")

    def __start_timer(self, name):
        self.__time_starts[name] = time.perf_counter()

    def __stop_timer(self, name):
        dt = time.perf_counter() - self.__time_starts[name]
        self.__time_elapsed[name] = self.__time_elapsed.get(name, 0.0) + dt

    def __get_list_k_plus_q(self, q):
        """return list of k-indices after k+q operation  (list of size nk_all)"""

        def k_plus_q(_k, _q):
            return [(_k[l] + _q[l]) % self.div[l] for l in range(3)]

        def k2i(_k):
            return ( _k[0] * self.div[1] + _k[1] ) * self.div[2] + _k[2]

        return [ k2i(k_plus_q(k, q)) for k in product(list(range(self.div[0])), list(range(self.div[1])), list(range(self.div[2]))) ]

    def __shift_by_q(self, g_k, q=(0, 0, 0)):
        """return G(k+q)  (np.ndarray of size nk_all)"""

        assert g_k.shape == (self.nk_all,)
        index_kq = self.__get_list_k_plus_q(q)
        g_kq = np.array([g_k[index_kq[j]] for j in range(self.nk_all)])
        return g_kq

    def __calc_chi0_fixed_freq(self, g1_k, g2_k, q):
        """
        Calculate chi_0(q) = - ave_k { g1(k) * g2(k+q) }

        g1_k[ik] : ik in FBZ
        g2_k[ik] : ik in FBZ
        """

        g2_k_q = self.__shift_by_q(g2_k, q)  # g2(k+q)

        chi0 = -np.dot(g1_k, g2_k_q) / self.nk_all # should cast explicitly?
        return chi0

    def calc_chi0_sum(self, g1_wk, g2_wk, wb, q, _n_wf):
        """
        Calculate chi_0(iwb,q; iwf) = -ave_k { g1(iwf, k) * g2(iwf+iwb, k+q) }

        return: np.array of size [2*_n_wf]
        g1_wk[iw][ik] : iw for positive and negative, ik in FBZ
        g2_wk[iw][ik] : the same as above
        wb : bosonic Matsubara frequency (>=0)
        q[kx,ky,kz] : momentum transfer
        _n_wf: cutoff for iwf
        """

        chi0 = np.array([ self.__calc_chi0_fixed_freq(g1_wk[wf], g2_wk[wf+wb], q) for wf in range(2*_n_wf) ])
        return chi0

    def __fft_k2r(self, gwk):
        gwk_xyz = np.reshape( gwk, (gwk.shape[0], self.div[0], self.div[1], self.div[2]) )
        gwr_xyz = np.fft.ifftn(gwk_xyz, axes=(1,2,3))
        gwr = np.reshape( gwr_xyz, (gwk.shape[0], self.nk_all,) )
        return gwr

    def __fft_r2k(self, gwr):
        gwr_xyz = np.reshape( gwr, (gwr.shape[0], self.div[0], self.div[1], self.div[2]) )
        gwk_xyz = np.fft.fftn(gwr_xyz, axes=(1,2,3))
        gwk = np.reshape( gwk_xyz, (gwr.shape[0], self.nk_all,) )
        return gwk

    def __calc_minus_r(self):
        nx = self.div[0]
        ny = self.div[1]
        nz = self.div[2]
        mr = np.empty((nx,ny,nz), dtype=int)
        for r,(x,y,z) in enumerate(product(list(range(nx)),list(range(ny)),list(range(nz)))):
            mr[(nx-x)%nx, (ny-y)%ny, (nz-z)%nz] = r
        return mr.reshape((nx*ny*nz))

    def calc_chi0_fft(self, g1_wk, g2_wk, wb, _n_wf):
        """return chi0_wq of size [2*_n_wf][nk_all]"""

        g1_wr = self.__fft_k2r(g1_wk)
        g2_wr = self.__fft_k2r(g2_wk)

        # -g1(w, r) * g2(w+w', -r)
        g2_wr_inversion = g2_wr[:, self.minus_r]  # g2(w, -r)
        chi0_wr = -g1_wr[0:2*_n_wf, :] * g2_wr_inversion[wb:wb+2*_n_wf, :]
        chi0_wq = self.__fft_r2k(chi0_wr)
        return chi0_wq

    def __get_G_iw_k(self, bname, icrsh1, icrsh2, in1, in2, iw_cutoff=(None,None), temp_file_ar=None):
        """
        return G(iw, k) for a given (bname, icrsh1, icrsh2, in1, in2)

        If iw_cutoff=(wmin,wmax) is given, the range of Matsubara freq is reduced to n=[wmin:wmax)
        """

        wmin, wmax = get_freq_index(n_iw=self.Sigma_imp_iw[0][bname].data.shape[0], iw_cutoff=iw_cutoff)

        # return g(iw,k)
        if temp_file_ar is None:
            # Gk_corr : [bname][icrsh, icrsh, in1, in2, iw, ik]
            return self.Gk_corr[bname][icrsh1, icrsh2, in1, in2, wmin:wmax, :]
        else:
            ds_name = 'crsh{}_crsh{}_{}_in{}_in{}'.format(icrsh1, icrsh2, bname, in1, in2)
            return temp_file_ar[ds_name][()][wmin:wmax, :]

    # TODO: offdiagonal (icrsh1, icrsh2)
    #       extract_G_loc needs to be extended
    def __get_G_loc_iw(self, bname, icrsh1, icrsh2, in1, in2, iw_cutoff=(None,None)):
        """
        return G_loc(iw) for a given (icrsh, bname, in1, in2)

        If iw_cutoff=(wmin,wmax) is given, the range of Matsubara freq is reduced to n=[wmin:wmax)
        """

        wmin, wmax = get_freq_index(n_iw=self.G_loc_corr[icrsh1][icrsh2][bname].data.shape[0], iw_cutoff=iw_cutoff)

        # g(iw)
        g_iw = np.array( self.G_loc_corr[icrsh1][icrsh2][bname].data[wmin:wmax, in1, in2] )
        return g_iw

    def __reduce_Sigma_iw(self, n_cutoff):
        """Reduce number of Matsubara frequencies in self.Sigma_imp_iw"""

        for icrsh in range(self.n_corr_shells):
            self.Sigma_imp_iw[icrsh] = gf_reduce_iw(self.Sigma_imp_iw[icrsh], n_cutoff)

        # delete cache in lattice_gf()
        if hasattr(self, "G_latt_iw"):
            del self.G_latt_iw

    def save_X0q_for_bse(
        self,
        list_wb,
        n_wf_cutoff,
        h5_file="dmft_bse.h5",
        groupname="",
        h5_mode="a",
        h5_compression="gzip",
        h5_compression_opts=None,
        with_Sigma=True,
        algo="fft",
        q_dict=None,
        flag_save_X0_loc=True,
        flag_save_chi0_loc=True,
        qpoints_saved="quadrant",
        nonlocal_order_parameter=False,
        temp_file=None,
        del_hopping=True,
    ):
        """
        list_wb: (list) list of bosonic frequencies (int) computed
        n_wf_cutoff: (int) cutoff for fermionic frequencies stored (2*n_wf_cutoff points are saved)
        del_hopping: (bool) If this is True, hopping is deleted from the memory. Subsequent call of other methods may fail.

        See also _save_X0_info, _save_X0_data
        """
        # TODO: difference between n_wf_store and n_wf_reduce?

        assert isinstance(list_wb, list)
        assert isinstance(n_wf_cutoff, int)

        wt = WallTime(str_time="\n Time:")

        # compute G_loc[icrsh1][icrsh2] using SK_ibz (defined with IBZ data)
        if mpi.is_master_node():
            print("\nCompute G_loc_corr")

        # first set necessary quantities to SK_ibz
        self.SK_ibz.set_mu(self.chemical_potential)
        self.SK_ibz.set_dc(self.dc_imp, self.dc_energ)

        # Sigma_imp_iw is defined for correlated shells
        assert len(self.Sigma_imp_iw) == self.n_corr_shells
        # Convert Sigma_imp_iw to inequivalent shells and pass it to set_Sigma
        Sigma_imp_iw_sh = [self.Sigma_imp_iw[self.inequiv_to_corr[ish]] for ish in range(self.n_inequiv_shells)]
        self.SK_ibz.set_Sigma(Sigma_imp_iw_sh)

        # then compute G_loc
        self.G_loc_corr = self.SK_ibz.extract_G_loc_offdiagonal(with_Sigma=with_Sigma)

        # TODO: check if Sigma_iw and G_loc_corr have the same structure

        # Reduce number of Matsubara frequecies
        # if n_wf_reduce is not None:
        if True:
            mpi.barrier()
            # self.__reduce_Sigma_iw(n_cutoff=n_wf_reduce)
            self.__reduce_Sigma_iw(n_cutoff=n_wf_cutoff+max(list_wb))

        if mpi.is_master_node():
            print("n_k =", self.n_k)
            def str_number_of_iw(_g):
                return ("%s"%_g.mesh).split(',')[0]
            print(str_number_of_iw(self.G_loc_corr[0][0]), " (G_loc)")  # print only # of iw
            print(str_number_of_iw(self.Sigma_imp_iw[0]), " (Sigma)")  # print only # of iw
            # print(("%s"%self.Sigma_imp_iw[0].mesh).split(',')[0], " (Sigma)")  # print only # of iw
            wt.print_time("G_loc_corr")

        # TODO: check if G_loc are consistent
        # self.G_loc_corr
        # self.extract_G_loc(with_Sigma=with_Sigma)

        # Compute Gk[ik][icrsh1][icrsh2]
        if mpi.is_master_node():
            print("\nCompute Gk_corr")
        if temp_file is None:
            self.Gk_corr = self.extract_Gk_crsh(with_Sigma=with_Sigma)
        else:
            self.extract_Gk_crsh(with_Sigma=with_Sigma, temp_file=temp_file)
        if mpi.is_master_node():
            wt.print_time("Gk_corr")

        # workaround to avoid memory error
        if mpi.is_master_node():
            print(f"\nMemory size of hopping:")
            print(f"  {self.hopping.nbytes:,} Bytes/process")
            print(f"  {self.hopping.nbytes*mpi.size:,} Bytes (total)")
            if del_hopping:
                print(f"del_hopping=True: hopping is deleted to save memory.")
            else:
                print(f"del_hopping=False: If this is True, hopping is deleted to save memory.")
        if del_hopping:
            del self.hopping

        # Compute X0(q) and save
        # self.init_for_X0(with_Sigma, n_wf_reduce)
        if mpi.is_master_node():
            print("\nCompute X0q")

            self.BS = h5BSE(
                h5_file,
                groupname,
                compression=h5_compression,
                compression_opts=h5_compression_opts,
            )
            self.BS.open(h5_mode)

        self._save_X0_info(nonlocal_order_parameter=nonlocal_order_parameter)
        for wb in list_wb:
            assert isinstance(wb, int)
            self._save_X0_data(wb, n_wf_cutoff, algo, q_dict, flag_save_X0_loc, qpoints_saved, temp_file=temp_file)
            if flag_save_chi0_loc:
                self._save_chi0_loc(wb)

        if mpi.is_master_node():
            self.BS.close()
            wt.print_time("X0q")

    def _save_X0_info(self, nonlocal_order_parameter):

        gf_struct = self.gf_struct_solver[0]  # ish=0
        for ish in range(1, self.n_inequiv_shells):
            assert gf_struct == self.gf_struct_solver[ish], "gf_struct should be identical among all correlated shells"
        block_g = list(gf_struct.keys())
        inner_g = gf_struct[block_g[0]]
        for inner in list(gf_struct.values()):
            assert inner == inner_g, "The inner structure should be identical among all blocks"

        for bname, g in self.Sigma_imp_iw[0]:
            beta = g.mesh.beta
            # n_wf = g.mesh.last_index() - g.mesh.first_index() + 1
            # positive_only = g.mesh.positive_only()  # this works only in ver1.5
            # TODO: check positive_only==False
            positive_only=False  #************ temporary
            break

        # only_diagonal1=True --> only local order parameter
        #                         = exclude non-local (unconventional) order parameter
        only_diagonal = not nonlocal_order_parameter

        self.block2 = IndexPair2(list(range(self.n_corr_shells)), block_g, only_diagonal1=only_diagonal)
        self.inner2 = IndexPair(inner_g, convert_to_int=True)

        if mpi.is_master_node():
            print(" G1 block:", block_g)
            print(" G1 inner:", inner_g)
            print(" G2 block:", self.block2.namelist)
            print(" G2 inner:", self.inner2.namelist)

        assert positive_only == False

        # save info
        self.X0_info = {}
        self.X0_info['block_name'] = self.block2.namelist
        self.X0_info['inner_name'] = self.inner2.namelist
        self.X0_info['beta'] = beta

        if mpi.is_master_node():
            self.BS.save( key=('block_name'), data=self.block2.namelist )
            self.BS.save( key=('inner_name'), data=self.inner2.namelist )
            self.BS.save( key=('beta'), data=beta )

    def __list_qpoints_saved(self, qpoints_saved):
        """return a list of ((int)q, (string)label) to be saved"""
        def q2str(_q):
            # return "%04d" %j
            qx = _q // (self.div[1]*self.div[2])
            res = _q % (self.div[1]*self.div[2])
            qy = res // self.div[2]
            qz = res % self.div[2]
            return "%02d.%02d.%02d" %(qx,qy,qz)

        def flatten_q(q_tuple):
            assert isinstance(q_tuple, tuple)
            qx, qy, qz = q_tuple
            return qz + self.div[2] * (qy + self.div[1] * qx)

        if isinstance(qpoints_saved, list):
            qpoints = list(map(flatten_q, qpoints_saved))
            # Remove duplicate elements (use sort to make the list unique order)
            qpoints = sorted(set(qpoints))
            q_labels = [ (q, q2str(q)) for q in qpoints]
        elif qpoints_saved == 'fbz':  # all q
            q_labels = [ (q, q2str(q)) for q in range(self.nk_all) ]
        # elif qpoints_saved == 'ibz':  # q in IBZ
        #     q_labels = []
        #     for kx,ky,kz in self.kpoint_ibz:
        #         q = (kx * self.div[1] + ky ) * self.div[2] + kz
        #         # print q, kx, ky, kz
        #         q_labels.append( (q, q2str(q)) )
        elif qpoints_saved == 'quadrant':  # q in the first quadrant
            q_labels = []
            for kx,ky,kz in product(list(range(self.div[0]//2+1)), list(range(self.div[1]//2+1)), list(range(self.div[2]//2+1))):
                q = (kx * self.div[1] + ky ) * self.div[2] + kz
                q_labels.append( (q, q2str(q)) )
        else:
            raise ValueError("_save_X0_data: qpoints_saved = '%s' not recognized" %qpoints_saved)
        return q_labels

    def __save_X0_data_sum(self, wb, n_wf_store, q_dict, temp_file_ar=None):
        # MPI not supported yet
        if not mpi.is_master_node():
            return

        for q_label, q_vec in list(q_dict.items()):
            chi0 = {}  # {(block, block) : chi}
            for (b_l, (icrsh1, b1, icrsh2, b2)), (b_r, (icrsh3, b3, icrsh4, b4)) in product(enumerate(self.block2.pairlist), repeat=2):
                # b1==b3, b2==b4
                if b1 != b3 or b2 != b4:
                    continue
                # print("", self.block2.namelist[b_l], "", self.block2.namelist[b_r])

                chi0_ijw = np.empty((self.inner2.size, self.inner2.size, 2 * n_wf_store), dtype=complex)
                for (i, (i1, i2)), (j, (i3, i4)) in product(enumerate(self.inner2.pairlist), repeat=2):
                    # print(" ", i1, i2, i3, i4)
                    g_l = self.__get_G_iw_k(b1, icrsh3, icrsh1, i3, i1, iw_cutoff=(-n_wf_store, n_wf_store + wb), temp_file_ar=temp_file_ar)
                    g_r = self.__get_G_iw_k(b2, icrsh2, icrsh4, i2, i4, iw_cutoff=(-n_wf_store, n_wf_store + wb), temp_file_ar=temp_file_ar)
                    temp = self.calc_chi0_sum(g_l, g_r, wb, q_vec, n_wf_store)  # [wf]
                    chi0_ijw[i, j, :] = temp
                chi0[(b_l, b_r)] = chi0_ijw

            if mpi.is_master_node():
                self.BS.save(key=('X0_q', wb, q_label), data=chi0)

    def __save_X0_data_fft(self, wb, n_wf_store, flag_save_X0_loc, qpoints_saved, temp_file_ar=None):

        self.__start_timer("save_X0_data_fft")

        # q-points to be saved
        #   list of tuple ((int)q, (str)label)
        q_labels = self.__list_qpoints_saved(qpoints_saved)

        # 1D array of sets of inner indices
        #   list of ( (i, (i1, i2)), (j, (i3, i4)) )
        ij_array = list(product(enumerate(self.inner2.pairlist), repeat=2))
        ij_indices = mpi.slice_array(np.arange(len(ij_array)))
        ij_num = len(ij_indices)

        data_sliced = [
            np.empty((ij_num, 2 * n_wf_store), dtype=np.complex128)
            for q in range(self.nk_all)
        ]  # [q][ [i,j,[w]] ]
        qsum_sliced = np.empty((ij_num, 2*n_wf_store), dtype=np.complex128)
        for (b_l, (icrsh1, b1, icrsh2, b2)), (b_r, (icrsh3, b3, icrsh4, b4)) in product(
            enumerate(self.block2.pairlist), repeat=2
        ):
            # b1==b3, b2==b4
            if b1 != b3 or b2 != b4:
                continue

            for idx, m in enumerate(ij_indices):
                (i, (i1, i2)), (j, (i3, i4)) = ij_array[m]
                # print(" rank%-3d:" % mpi.rank, i, j, i1, i2, i3, i4)
                g_l = self.__get_G_iw_k(b1, icrsh3, icrsh1, i3, i1, iw_cutoff=(-n_wf_store, n_wf_store + wb), temp_file_ar=temp_file_ar)
                g_r = self.__get_G_iw_k(b2, icrsh2, icrsh4, i2, i4, iw_cutoff=(-n_wf_store, n_wf_store + wb), temp_file_ar=temp_file_ar)
                chi0_wq = self.calc_chi0_fft(g_l, g_r, wb, n_wf_store)

                # store only part of q component
                for q, label in q_labels:
                    data_sliced[q][idx,:] = chi0_wq[:, q]
                qsum_sliced[idx,:] = np.average(chi0_wq, axis=1)  # average over q

                # chi0_wq must be deleted to release memory
                del g_l, g_r, chi0_wq

            # gather chi0_q data
            for q, label in q_labels:
                send_buf = data_sliced[q].reshape(-1)
                data_gathered = gatherv(send_buf)

                # save data to HDF5
                if mpi.is_master_node():
                    chi0_ijw = data_gathered.reshape((self.inner2.size, self.inner2.size, 2 * n_wf_store))

                    self.__start_timer("save")
                    self.BS.save(key=('X0_q', wb, label), data={(b_l, b_r): chi0_ijw})
                    self.__stop_timer("save")
                # mpi.barrier()
                # del data_sliced[q]

            # gather chi0_loc
            send_buf = qsum_sliced.reshape(-1)
            qsum_gathered = gatherv(send_buf)
            if mpi.is_master_node():
                chi0_qsum = qsum_gathered.reshape((self.inner2.size, self.inner2.size, 2 * n_wf_store))

                # save data to HDF5
                if flag_save_X0_loc:
                    self.__start_timer("save")
                    self.BS.save(key=('X0_loc', wb), data={(b_l, b_r): chi0_qsum})
                    self.__stop_timer("save")
        self.__stop_timer("save_X0_data_fft")

    def _save_X0_data(self, wb, n_wf_store, algo='fft', q_dict=None, flag_save_X0_loc=True, qpoints_saved='quadrant',
                     temp_file=None):
        """
        wb: (int) bosonic frequency
        algo: (string) how to compute X0_q: 'fft' for FFT, and 'sum' for k-sum
        q_dict: (dict) q vectors to be computed when algo='sum'. Ignored when algo='fft'
            key = q_label: (string) postfix of subgroup in hdf5.
            val = q_vec: (int, int, int) specific momentum to compute.
        n_wf_store: (int) cutoff for fermionic frequencies stored (2*n_wf_store points are saved)
        qpoints_saved: 'fbz', 'ibz', 'quandant'
        """

        if temp_file is None:
            temp_file_ar = None
        else:
            temp_file_ar = HDFArchive(os.path.abspath(temp_file), 'r')

        assert algo.lower() in ['fft', 'sum']
        if algo.lower() == 'sum':
            assert isinstance(q_dict, dict)
            if mpi.is_master_node():
                for q_label, q_vec in list(q_dict.items()):
                    print("# w=", wb, ", q=", [str(q_vec[l]) + "/" + str(self.div[l]) for l in range(3)], " (%s)" % q_label)
            self.__save_X0_data_sum(wb, n_wf_store, q_dict, temp_file_ar=temp_file_ar)

        else:  # fft
            if mpi.is_master_node():
                print("# w=", wb, ", q=all (using FFT)")
            self.__save_X0_data_fft(wb, n_wf_store, flag_save_X0_loc, qpoints_saved, temp_file_ar=temp_file_ar)

        mpi.barrier()

        if not temp_file_ar is None:
            del temp_file_ar

    def save_X_loc_dummy(self, wb, n_wf_store):
        """Save dummy data (zero matrix) for X_loc to enable running BSE solver without actually computing X_loc"""
        if not mpi.is_master_node():
            return

        inner_chi_name = self.X0_info['inner_name']
        X_loc_dummy = np.zeros( (len(inner_chi_name),len(inner_chi_name),2*n_wf_store,2*n_wf_store), dtype=complex )
        self.BS.save( key=('X_loc', wb), data={(0,0) : X_loc_dummy} )

    def _save_chi0_loc(self, wb):
        """Compute q-averaged susceptibility by taking summation over fermionic Matsubara frequency"""
        if not mpi.is_master_node():
            return

        beta = self.X0_info['beta']

        suscep = {}  # {(block, block) : chi}
        for (b_l, (icrsh1, b1, icrsh2, b2)), (b_r, (icrsh3, b3, icrsh4, b4)) in product(enumerate(self.block2.pairlist), repeat=2):
            # b1==b3, b2==b4
            if b1 != b3 or b2 != b4:
                continue
            # print("", self.block2.namelist[b_l], "", self.block2.namelist[b_r])

            suscep_ij = np.empty( (self.inner2.size, self.inner2.size), dtype=complex )
            for (i, (i1, i2)), (j, (i3, i4)) in product(enumerate(self.inner2.pairlist), repeat=2):
                # print(" ", i1, i2, i3, i4)
                g_l = self.__get_G_loc_iw(b1, icrsh3, icrsh1, i3, i1)
                g_r = self.__get_G_loc_iw(b2, icrsh2, icrsh4, i2, i4)
                wmax = g_l.shape[0]
                if wb>=0:
                    suscep_ij[i,j] = (-1./beta) * np.dot(g_l[:wmax-wb], g_r[wb:])  # sum_wf g(wf) * g(wf+wb)
                else:
                    suscep_ij[i,j] = (-1./beta) * np.dot(g_l[-wb:], g_r[:wmax+wb])
            suscep[(b_l, b_r)] = suscep_ij

        self.BS.save( key=('chi0_loc_in', wb), data=suscep )

    def extract_Gk_crsh(self, mu=None, iw_or_w='iw', with_Sigma=True, with_dc=True, broadening=None, temp_file=None):
        r"""
        [ORIGINAL: function 'extract_G_loc' in sumk_dft.py]
        Extracts the k-resolved downfolded Green function.

        Parameters
        ----------
        mu : real, optional
             Input chemical potential. If not provided the value of self.chemical_potential is used as mu.
        with_Sigma : boolean, optional
                     If True then the local GF is calculated with the self-energy self.Sigma_imp.
        with_dc : boolean, optional
                  If True then the double-counting correction is subtracted from the self-energy in calculating the GF.
        broadening : float, optional
                     Imaginary shift for the axis along which the real-axis GF is calculated.
                     If not provided, broadening will be set to double of the distance between mesh points in 'mesh'.
                     Only relevant for real-frequency GF.

        Returns
        -------
        G_k : dict of numpy.ndarray, [bname][icrsh, icrsh, in1, in2, iw, ik]
              k-dependent Green's functions
              rotated into the corresponding local frames.

        """

        # icrsh=0
        mesh = self.Sigma_imp_iw[0].mesh
        beta = mesh.beta
        bnames = [bname for bname, g in self.Sigma_imp_iw[0]]
        bname0 = bnames[0]
        inner = self.Sigma_imp_iw[0][bname0].indices
        inner_size = get_block_size(self.Sigma_imp_iw[0][bname0])
        n_iw = self.Sigma_imp_iw[0][bname0].data.shape[0]

        # Sigma_imp_iw is defined for correlated shells
        assert len(self.Sigma_imp_iw) == self.n_corr_shells

        # All correlated shells have the same block/inner structure
        for icrsh in range(self.n_corr_shells):
            assert mesh == self.Sigma_imp_iw[icrsh].mesh

        if mu is None:
            mu = self.chemical_potential

        # G_loc is rotated to the local coordinate system:
        # if self.use_rotations:
        #     for r in self.rot_mat:
        #         if not np.allclose(np.identity(r.shape[0]), r):
        #             raise Exception("Not implemented: rot_mat should an identity matrix!")

        def rotate_gf(gf, icrsh1, icrsh2):
            if self.use_rotations:
                gf.from_L_G_R(self.rot_mat[icrsh1].conjugate().transpose(), gf, self.rot_mat[icrsh2])

        if iw_or_w != "iw":
            raise Exception("Only iw_or_w=='iw' implemented")

        # Gk_corr : [bname][icrsh, icrsh, in1, in2, iw, ik]
        gk_shape = (self.n_corr_shells, self.n_corr_shells, inner_size, inner_size, n_iw, self.n_k)
        def init_Gk_corr():
            return {bname: np.zeros(gk_shape, dtype=complex) for bname in bnames}

        if mpi.is_master_node():
            z = np.zeros(1, dtype=complex)
            # print(z.nbytes)
            gk_nbytes = np.prod(gk_shape) * len(bnames) * z.nbytes
            print(f" Memory size of Gk_corr: {gk_nbytes:,} Bytes")
            print(f" | for each process if use_temp_files=False")
            print(f" | for all processes if True (each process store only a part of Gk_corr)")

        # g_template = BlockGf(name_block_generator=[(block, GfImFreq(indices=inner, mesh=mesh)) for block in bnames], make_copies=False)

        g_temp = GfImFreq(indices=inner, mesh=mesh)

        how_to_parallelize = 2
        # *****************************************************************
        # CALC G_k[ik][icrsh]
        # 0. No MPI
        # *****************************************************************
        if how_to_parallelize==0:
            assert temp_file is None
            Gk_corr = init_Gk_corr()
            for ik in range(self.n_k):
                G_latt = self.lattice_gf(ik=ik, mu=mu, iw_or_w=iw_or_w, with_Sigma=with_Sigma, with_dc=with_dc, beta=beta)

                for bname, gf in G_latt:
                    for icrsh1, icrsh2 in product(np.arange(self.n_corr_shells), repeat=2):
                        g_temp << self.downfold_offdiagonal(ik, icrsh1, icrsh2, bname, gf, g_temp)
                        rotate_gf(g_temp, icrsh1, icrsh2)
                        # [iw, in1, in2]
                        g_wij = g_temp.data
                        # [in1, in2, iw]
                        Gk_corr[bname][icrsh1, icrsh2, :, :, :, ik] = g_wij.transpose((1, 2, 0))  # copy
            # return Gk_corr

        # *****************************************************************
        # 2. MPI. use Allgather to treat large-scale data
        # *****************************************************************
        elif how_to_parallelize==2:

            ikarray = np.arange(self.n_k)
            ikarray_sliced = mpi.slice_array(ikarray)

            # Test if gathering sliced arrays reproduces the original array
            def test_gather(array):
                """test for slice_array and MPI.COMM_WORLD.allgather"""
                array_gather = MPI.COMM_WORLD.allgather(ikarray_sliced)
                array_gather = np.concatenate(array_gather, axis=0)
                return np.allclose(array, array_gather)
            assert test_gather(ikarray)

            # [ik, icrsh, icrsh, in1, in2, iw]
            gk_sliced_shape = (ikarray_sliced.size,) + gk_shape[:-1]
            gk_sliced = {bname: np.zeros(gk_sliced_shape, dtype=complex) for bname in bnames}

            for jk, ik in enumerate(ikarray_sliced):
                G_latt = self.lattice_gf(ik=ik, mu=mu, iw_or_w=iw_or_w, with_Sigma=with_Sigma, with_dc=with_dc, beta=beta)

                for bname, gf in G_latt:
                    for icrsh1, icrsh2 in product(np.arange(self.n_corr_shells), repeat=2):
                        g_temp << self.downfold_offdiagonal(ik, icrsh1, icrsh2, bname, gf, g_temp)
                        rotate_gf(g_temp, icrsh1, icrsh2)
                        # [iw, in1, in2]
                        g_wij = g_temp.data
                        # [in1, in2, iw]
                        gk_sliced[bname][jk, icrsh1, icrsh2, :, :, :] = g_wij.transpose((1, 2, 0))  # copy

            if temp_file is None:
                Gk_corr = init_Gk_corr()
                for bname in bnames:
                    gk_gather = MPI.COMM_WORLD.allgather(gk_sliced[bname])
                    gk_gather = np.concatenate(gk_gather, axis=0)

                    # [ik, icrsh, icrsh, in1, in2, iw]
                    assert gk_gather.shape == (self.n_k,) + gk_shape[:-1]

                    # [icrsh, icrsh, in1, in2, iw, ik]
                    Gk_corr[bname][...] = gk_gather.transpose((1, 2, 3, 4, 5, 0))  # copy

                del gk_gather
                # return Gk_corr
            else:
                if mpi.rank == 0:
                    ar = HDFArchive(os.path.abspath(temp_file), 'w')
                for bname in bnames:
                    for icrsh1, icrsh2 in product(np.arange(self.n_corr_shells), repeat=2):
                        for in1, in2 in product(np.arange(inner_size), repeat=2):
                            gk_gather = MPI.COMM_WORLD.gather(gk_sliced[bname][:, icrsh1, icrsh2, in1, in2, :], root=0)
                            if mpi.rank == 0:
                                gk_gather = np.concatenate(gk_gather, axis=0)
                                ds_name = 'crsh{}_crsh{}_{}_in{}_in{}'.format(icrsh1, icrsh2, bname, in1, in2)
                                # save values of the array as (iwn, k)
                                ar[ds_name] = gk_gather.transpose()
                if mpi.rank == 0:
                    del ar
                    # TODO print file size of temp_file
                mpi.barrier()
                Gk_corr = None

            del gk_sliced

        return Gk_corr
