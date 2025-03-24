#!/usr/bin/env python3

import numpy as np
from itertools import product
import os
import math
import sys
import argparse
import configparser

from chiq.h5bse import h5BSE
from chiq.matrix_dict import MatrixDict
from chiq.g2scl_core import G2SCL_core

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import cm
# from mpl_toolkits.mplot3d import Axes3D


def plot_2d_array(array, filename, ylabel='', legend=None):
    assert array.ndim == 2
    # plt.plot(array[:,0], array[:,1], marker='.')
    # plt.plot(array[:,0], array[:,1], 'o')
    for y in range(1, array.shape[1]):
        plt.plot(array[:,0], array[:,y], 'o', label=legend[y-1] if legend is not None else '')
    plt.legend()
    plt.ylabel(ylabel)
    plt.axhline(0) # y=0
    plt.savefig(filename)
    print(" '%s'" % filename)
    plt.close()


def save_array(array, filename, flag_index=True):
    assert array.ndim == 1
    np.savetxt(filename, np.array(list(enumerate(array))), fmt=['%d', '%.18e'])
    print(" '%s'" % filename)


def symmetrize(x_loc):
    # TODO
    raise NotImplementedError
    # return x_loc, 0


def get_shape(_x_dict):
    x_keys = list(_x_dict.keys())
    x_shape = _x_dict[x_keys[0]].shape
    # check if all components in dict have the same shape
    for key, array in list(_x_dict.items()):
        # print array.dtype
        # assert isinstance(array.dtype, np.complex128)
        assert array.shape == x_shape
    return x_keys, x_shape


class G2SCL:
    # def __init__(self, X_dict, chi_dict, algo, beta, flag_plot=True, flag_3dplot=True, flag_save_3ddata=False, dir_plot='scl', flag_s0_approx=False, handle_exception='zero', **kwargs):
    def __init__(self, beta, flag_plot=True, flag_3dplot=True, flag_save_3ddata=False, dir_plot='scl', flag_s0_approx=False, handle_exception='zero', **kwargs):

        """

        beta : [float] only for plot

        """

        self.beta = beta
        self.flag_plot = flag_plot
        self.flag_3dplot = flag_3dplot and flag_plot  # flag_3dplot becomes False if flag_plot is False
        self.flag_save_3ddata = flag_save_3ddata
        self.dir_plot = dir_plot
        self.handle_exception = handle_exception
        self.flag_s0_approx = flag_s0_approx

        if not os.path.exists(dir_plot):
            os.mkdir(dir_plot)

        # to be set in set_X()
        self.algo = None
        self.x_dict_reshape = None
        self.Nw = None
        self.Nw_reshape = None
        self.opt = {}  # optional variables dependent on algo

        # to be set in set_eigenmodes()
        self.M = None
        self.n_mode = None
        self._md = None
        self.md_keys = None
        self._flag_eigenmodes = False

        # to be set in chi_loc()
        self.chi_loc = None

        # to be set in compute()
        self.phi_mode = None

    def set_x(self, algo, X_dict, **kwargs):
        """

        Parameters
        ----------
        algo : [str]
        X_dict   : [dict] dictionary block matrix of two-particle GF, X;
                   Each block contain numpy.ndarray whose shape depends on algo.
        """
        print("\nSet X_loc")
        self._check_if_set_eigenmodes_called()

        self.algo = algo

        # check data structure of X_dict
        def get_x_shape(_x_dict):
            _x_keys, _x_shape = get_shape(_x_dict)
            assert _x_keys == self.md_keys
            assert _x_shape[0] == _x_shape[1] == self.M
            return _x_shape

        if algo == 'svd':
            # [M, M, Nw, Nw]
            x_shape = get_x_shape(X_dict)
            assert len(x_shape) == 4
            assert x_shape[2] == x_shape[3]
            self.Nw = x_shape[2]
            # [M, M, Nw*Nw]
            self.Nw_reshape = self.Nw*self.Nw
            self.x_dict_reshape = {key: np.reshape(array, (self.M, self.M, self.Nw_reshape))
                                   for key, array in list(X_dict.items())}
        elif algo == '1dhalf':
            # [M, M, Nw/2]
            x_shape = get_x_shape(X_dict)
            assert len(x_shape) == 3
            self.Nw = x_shape[2]*2
            # self.n = kwargs['n']
            # self.z = kwargs['z']
            self.opt['n'] = kwargs['n']
            self.opt['z'] = kwargs['z']
            # [M, M, Nw/2]
            self.Nw_reshape = self.Nw // 2
            self.x_dict_reshape = X_dict
        elif algo == '1dfull':
            # [M, M, Nw]
            x_shape = get_x_shape(X_dict)
            assert len(x_shape) == 3
            self.Nw = x_shape[2]
            # self.n = kwargs['n']
            self.opt['n'] = kwargs['n']
            # [M, M, Nw]
            self.Nw_reshape = self.Nw
            self.x_dict_reshape = X_dict
        elif algo == 'anti_diag':
            # [M, M, Nw/2]
            x_shape = get_x_shape(X_dict)
            assert len(x_shape) == 3
            self.Nw = x_shape[2]*2
            # [M, M, Nw/2]
            self.Nw_reshape = self.Nw // 2
            self.x_dict_reshape = X_dict
        elif algo == '2pole':
            # self.x_dict_reshape = None
            self.Nw = kwargs['iw_cutoff']*2
            # self.E_plus = kwargs['E_plus']
            # self.E_minus = kwargs['E_minus']
            self.opt['E_plus'] = kwargs['E_plus']
            self.opt['E_minus'] = kwargs['E_minus']
        else:
            raise ValueError("algo=%s not found" % algo)

    def _check_if_set_eigenmodes_called(self):
        if self._flag_eigenmodes == False:
            print("set_eigenmodes() should be called in advance.")
            exit(1)

    def _print_eigenmodes(self, eigens, basename):
        save_array(eigens, '%s/%s.dat' %(self.dir_plot, basename))
        if flag_plot:
            plot_2d_array(np.array(list(enumerate(eigens))), '%s/%s.pdf' %(self.dir_plot, basename))
        print("mode  eigenvalue")
        for mode, x in enumerate(eigens):
            print("", mode, "", x)
            # assert x.imag < 1e-12  # NOTE: may not be fulfilled in CTQMC
            # TODO
            # assert x.real >= 0

    def set_eigenmodes(self, _matrix_dict):
        """
        Set eigenmodes by diagonalizing _matrix_dict.

        Parameters
        ----------
        _matrix_dict: [dict]

        """
        print("\nSet eigenmodes")
        print("  CHECK DEGENERACY OF EIGENVALUES")
        self._flag_eigenmodes = True

        # check the structure of _matrix_dict
        self.md_keys, md_shape = get_shape(_matrix_dict)
        print("structure of matrix_dict:")
        print("  | keys  =", self.md_keys)
        print("  | shape =", md_shape)
        assert md_shape[0] == md_shape[1]
        self.M = md_shape[0]

        self._md = MatrixDict(_matrix_dict, verbose=False)
        self.n_mode = self._md.get_eigen().shape[0]
        print("n_mode =", self.n_mode)
        if True:
            eigen = self._md.get_eigen()
            self._print_eigenmodes(eigen, 'eigenmodes')
            self._md.print_eigenvectors(print_as_matrix=True)

    def set_chiloc(self, chi_loc):
        print("\nSet chi_loc")
        self._check_if_set_eigenmodes_called()

        self.chiloc_eigen, err = self._md.compute_eigen(chi_loc)

        # TODO: check if chiloc_eigen is real

        print("***Error in orbital diagonalization = %.2e" % err)
        if True:
            eigen = np.real(self.chiloc_eigen)
            self._print_eigenmodes(eigen, 'chi_loc_diag')

    def compute(self):
        print("\nCompute")

        # diagonalize
        # [M, Nw, Nw]
        if self.x_dict_reshape is not None:
            x_mode = np.zeros((self.n_mode, self.Nw_reshape), dtype=np.complex128)
            error_orbital_diag = 0
            for wf in range(self.Nw_reshape):
                x = {key: array[:, :, wf] for key, array in list(self.x_dict_reshape.items())}
                x_mode[:, wf], error = self._md.compute_eigen(x)
                error_orbital_diag = max(error_orbital_diag, error)
            print("***Error in orbital diagonalization = %.2e" % error_orbital_diag)
            print("x_mode.shape =", x_mode.shape)
            assert self.n_mode == x_mode.shape[0]

        # compute Phi_L in each mode
        # [M, Nw]
        self.phi_mode = np.zeros((self.n_mode, self.Nw), dtype=np.complex128)
        print("phi_mode.shape =", self.phi_mode.shape)
        s0_list = []
        for mode in range(self.n_mode):
            g2scl = G2SCL_core(self.beta, handle_exception=handle_exception)
            if self.algo == 'svd':
                # reshape [Nw*Nw] to [Nw,Nw]
                x_mode_reshape = np.reshape(x_mode[mode], (self.Nw, self.Nw))
                g2scl.set_exact(x_mode_reshape)
            elif self.algo == '1dhalf':
                z_mode, _ = self._md.compute_eigen(self.opt['z'])  # discard error
                g2scl.set_1dhalf(x_mode[mode], self.opt['n'], z_mode[mode])
            elif self.algo == '1dfull':
                g2scl.set_1dfull(x_mode[mode], self.opt['n'])
            elif self.algo == 'anti_diag':
                g2scl.set_anti_diag(x_mode[mode])
            elif self.algo == '2pole':
                g2scl.set_2pole(self.opt['E_plus'], self.opt['E_minus'], self.Nw)
            # print u0.shape

            # g2scl.set_chiloc(self._md.get_eigen()[mode])
            g2scl.set_chiloc(self.chiloc_eigen[mode])
            # phi_mode[mode,:] = G2D.get_Phi_L()
            self.phi_mode[mode, :] = g2scl.get_Phi_L(use_chiloc=flag_s0_approx)
            s0_list.append([mode, g2scl.get_s0(), g2scl.get_s0(use_chiloc=True)])

            print("mode =", mode, ", s0 =", g2scl.get_s0())
            print("               ", g2scl.get_s0(use_chiloc=True), "(from chi_loc)")

            subdir_plot = self.dir_plot + '/%02d' % mode
            if flag_plot:
                if not os.path.exists(subdir_plot):
                    os.mkdir(subdir_plot)
                g2scl.save_sv('%s/sv.dat' % subdir_plot)
                g2scl.plot_sv('%s/sv.pdf' % subdir_plot)
                g2scl.plot_u0('%s/u0.pdf' % subdir_plot, data_out='%s/u0.dat' % subdir_plot)
            if flag_3dplot:
                g2scl.plot_x('%s/X.pdf' % subdir_plot, use_chiloc=flag_s0_approx, flag_data_out=flag_save_3ddata)
            del g2scl
        # print phi_mode.shape
        # np.savetxt('%s/s0.dat' %dir_plot, s0_list, fmt=["%d", "%.18e"])
        np.savetxt('%s/s0.dat' % self.dir_plot, s0_list, fmt=["%d", "%.18e", "%.18e"])
        if flag_plot:
            plot_2d_array(np.array(s0_list), '%s/s0.pdf' % self.dir_plot, ylabel='$s_0$', legend=['s0', 'from chi_loc'])

    def save(self, h5_file='dmft_bse.h5', groupname=''):
        print("\nSave")

        if self.phi_mode is None:
            raise Exception("phi_mode has not been computed. Call compute() before save()")

        # transform Phi_L into the original spin-orbital representation
        # [M, M, Nw]
        phi_shape = (self.M, self.M, self.Nw)
        print("phi_dict[].shape =", phi_shape)
        phi_dict = {key: np.zeros(phi_shape, dtype=np.complex128) for key in self.md_keys}
        for wf in range(self.Nw):
            for key, data in list(self._md.transform_diag(self.phi_mode[:, wf]).items()):
                assert key in phi_dict
                phi_dict[key][:, :, wf] = data[:, :]

        BS = h5BSE(h5_file, groupname)
        BS.save(key=('Phi', 0), data=phi_dict)
        del BS


if __name__ == '__main__':
    P = argparse.ArgumentParser()
    P.add_argument('-f', default='dmft_bse.h5', help="h5 filename. default=dmft_bse.h5")
    P.add_argument('--algo', default='svd', help="comma-separated list including svd, 1dhalf, 1dfull, 2pole, and/or anti_diag. default=svd")
    P.add_argument('--eigenmode', default='chi_loc', choices=['chi_loc', 'mix'], help="Definition of eigenmodes.")
    P.add_argument('--weight_odd', type=float, default=0.1, help="Weight of odd-frequency components to define eigenmodes when --eigenmode=mix")
    P.add_argument('--noplot', action='store_false', help="No figure is generated")
    P.add_argument('--no3dplot', action='store_false', help="No 3D figure is generated (activated by --noplot)")
    P.add_argument('--save_3ddata', action='store_true', help="Save 3D data (no concurrent use with --no3dplot)")
    P.add_argument('--s0_approx', action='store_true', help="Estimate s0 from chi_loc approximately")
    P.add_argument('--input', default='scl.in', help="Input file name (for algo=2pole)")
    P.add_argument('--symmetrize', action='store_true', help="Symmetrize X_loc (activate this for QMC data)")
    P.add_argument('--handle_exception', default='zero', help="How to handle exception: zero (default), ignore, exit.")
    args = P.parse_args()
    print(args)

    h5_file = args.f
    flag_plot = args.noplot
    flag_3dplot = args.no3dplot and args.noplot  # True if both are True
    flag_save_3ddata = args.save_3ddata
    flag_s0_approx = args.s0_approx
    input_file = args.input
    flag_symmetrize = args.symmetrize
    handle_exception = args.handle_exception
    eigenmode = args.eigenmode
    weight_odd = args.weight_odd

    list_algo = args.algo.lower().split(",")
    print("list_algo =", list_algo)

    BS = h5BSE(h5_file)
    x_loc = BS.get(key=('X_loc', 0))
    beta = BS.get(key=('beta',))
    print("x_loc.keys =", list(x_loc.keys()))

    chi_loc = BS.get(key=('chi_loc', 0))
    # compute chi_loc from x_loc if data do not exist in the database
    if not isinstance(chi_loc, dict):
        chi_loc = {key: np.sum(array, axis=(2, 3))/beta for key, array in list(x_loc.items())}
    del BS

    if flag_symmetrize:
        x_loc, err = symmetrize(x_loc)
        print(" symmetrization error = %.3e" % err)
        chi_loc, err = symmetrize(chi_loc)
        print(" symmetrization error = %.3e" % err)

    # compute odd-frequency components in X_loc
    if True:
        x_nw = x_loc[list(x_loc.keys())[0]].shape[2]
        fac_odd = np.array([-1,]*(x_nw//2) + [1,]*(x_nw//2), dtype=float)  # odd-frequency factor
        # fac_asym = np.ones((x_nw/2,))
        # print(fac_odd)

        def odd_sum(array4):
            return np.einsum("n,ijnm,m->ij", fac_odd, array4, fac_odd) / beta

        chi_odd = {key: odd_sum(array) for key, array in list(x_loc.items())}

    # diagonalize = 'chi_loc'
    if eigenmode == 'chi_loc':
        chi_to_be_diagonalized = chi_loc
    elif eigenmode == 'mix':
        chi_loc_mix = {key: chi_loc[key] * (1.0-weight_odd) + chi_odd[key] * weight_odd
                       for key, array in list(chi_loc.items())}
        chi_to_be_diagonalized = chi_loc_mix

    # common
    w0 = x_loc[list(x_loc.keys())[0]].shape[3] // 2

    def omega_n(_n):
        return w0 + _n  # w_n

    def minus_omega_n(_n):
        return w0 - _n - 1  # -w_n

    params = {
        'beta': beta,
        'flag_plot': flag_plot,
        'flag_3dplot': flag_3dplot,
        'flag_save_3ddata': flag_save_3ddata,
        'flag_s0_approx': flag_s0_approx,
        'handle_exception': handle_exception,
    }
    print("params =", params)

    # svd
    if 'svd' in list_algo:
        print("\n=====================================\nalgo='svd'")
        # scl = G2SCL(x_loc, chi_loc, algo='svd', beta=beta, h5_file=h5_file, groupname='scl_svd', flag_plot=flag_plot, flag_3dplot=flag_3dplot, flag_save_3ddata=flag_save_3ddata, flag_s0_approx=flag_s0_approx, handle_exception=handle_exception)
        scl = G2SCL(dir_plot='scl', **params)
        scl.set_eigenmodes(chi_to_be_diagonalized)
        scl.set_x(algo='svd', X_dict=x_loc)
        scl.set_chiloc(chi_loc)
        scl.compute()
        scl.save(h5_file=h5_file, groupname='scl_svd')

    # 1dhalf
    if '1dhalf' in list_algo:
        print("\n=====================================\nalgo='1dhalf'")
        # X(iw, -iw_n)
        n = 0
        w_n = omega_n(n)
        mw_n = minus_omega_n(n)
        m = 1
        w_m = omega_n(m)
        mw_m = minus_omega_n(m)
        # print w0, wn, mwn
        x_1d = {key: array[:, :, w0:, mw_n] for key, array in list(x_loc.items())}
        z = {key: array[:, :, w_m, mw_n]+array[:, :, mw_m, mw_n] for key, array in list(x_loc.items())}

        # scl = G2SCL(x_1d, chi_loc, algo='1dhalf', beta=beta, h5_file=h5_file, groupname='scl_1dhalf', flag_plot=flag_plot, flag_3dplot=flag_3dplot, flag_save_3ddata=flag_save_3ddata, dir_plot='scl_1dhalf', flag_s0_approx=flag_s0_approx, handle_exception=handle_exception, n=n, z=z)
        scl = G2SCL(dir_plot='scl_1dhalf', **params)
        scl.set_eigenmodes(chi_to_be_diagonalized)
        scl.set_x(algo='1dhalf', X_dict=x_1d, n=n, z=z)
        scl.set_chiloc(chi_loc)
        scl.compute()
        scl.save(h5_file=h5_file, groupname='scl_1dhalf')

    # 1dfull
    if '1dfull' in list_algo:
        print("\n=====================================\nalgo='1dfull'")
        # X(iw, iw_n)
        n = 0
        x_1d = {key: array[:, :, :, omega_n(n)] for key, array in list(x_loc.items())}
        # scl = G2SCL(x_1d, chi_loc, algo='1dfull', beta=beta, h5_file=h5_file, groupname='scl_1dfull', flag_plot=flag_plot, flag_3dplot=flag_3dplot, flag_save_3ddata=flag_save_3ddata, dir_plot='scl_1dfull', flag_s0_approx=flag_s0_approx, handle_exception=handle_exception, n=n)
        scl = G2SCL(dir_plot='scl_1dfull', **params)
        scl.set_eigenmodes(chi_to_be_diagonalized)
        scl.set_x(algo='1dfull', X_dict=x_1d, n=n)
        scl.set_chiloc(chi_loc)
        scl.compute()
        scl.save(h5_file=h5_file, groupname='scl_1dfull')

    # anti_diag
    if 'anti_diag' in list_algo:
        print("\n=====================================\nalgo='anti_diag'")
        # X(iw,-iw) for w>0

        def func(array):
            return np.einsum("ijkk->ijk", array[:, :, w0:, w0-1::-1])

        x_anti_diag = {key: func(array) for key, array in list(x_loc.items())}
        # scl = G2SCL(x_anti_diag, chi_loc, algo='anti_diag', beta=beta, h5_file=h5_file, groupname='scl_antidiag', flag_plot=flag_plot, flag_3dplot=flag_3dplot, flag_save_3ddata=flag_save_3ddata, flag_s0_approx=flag_s0_approx, handle_exception=handle_exception, dir_plot='scl_antidiag')
        scl = G2SCL(dir_plot='scl_antidiag', **params)
        scl.set_eigenmodes(chi_to_be_diagonalized)
        scl.set_x(algo='anti_diag', X_dict=x_anti_diag)
        scl.set_chiloc(chi_loc)
        scl.compute()
        scl.save(h5_file=h5_file, groupname='scl_antidiag')

    if '2pole' in list_algo:
        print("\n=====================================\nalgo='2pole'")
        # read parameters in a file
        config = configparser.ConfigParser()
        config.read(input_file)
        E_plus_ = float(config.get('SCL', 'E_plus'))
        E_minus_ = float(config.get('SCL', 'E_minus'))
        iw_cutoff_ = int(config.get('SCL', 'iw_cutoff'))
        del config

        # scl = G2SCL(None, chi_loc, algo='2pole', beta=beta, h5_file=h5_file, groupname='scl_2pole', flag_plot=flag_plot, flag_3dplot=flag_3dplot, flag_save_3ddata=flag_save_3ddata, flag_s0_approx=True, handle_exception=handle_exception, dir_plot='scl_2pole', E_plus=E_plus_, E_minus=E_minus_, iw_cutoff=iw_cutoff_)
        scl = G2SCL(dir_plot='scl_2pole', **params)
        scl.set_eigenmodes(chi_to_be_diagonalized)
        scl.set_x(algo='2pole', X_dict=None, E_plus=E_plus_, E_minus=E_minus_, iw_cutoff=iw_cutoff_)
        scl.set_chiloc(chi_loc)
        scl.compute()
        scl.save(h5_file=h5_file, groupname='scl_2pole')

    if flag_plot:
        pass
