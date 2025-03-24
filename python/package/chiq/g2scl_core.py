

import numpy as np
import math
import os
from scipy.optimize import curve_fit, leastsq
import inspect

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D


def is_real(array):
    flag_real = True
    for z in array:
        if z.imag > 1e-10:  flag_real = False
    return flag_real


def phase(z):
    return z / abs(z)


def shift_phase(array, z):
    assert isinstance(array, np.ndarray)
    # fac = z / abs(z)  # phase factor
    return array*phase(z)


def matsubara_mesh(y_data, beta):
    """return Matsubara frequency"""
    w0 = y_data.shape[0] // 2
    return np.array([ (2*i+1)*math.pi/beta for i in range(-w0, w0) ])


def lorenzian(x, a, width):
    return a * width / (x**2 + width**2)


def fit_sym(x, y):
    print("  fit with lorenzian")
    try:
        params, cov = curve_fit(lorenzian, x, y)
        print("  | a =", params[0])
        print("  | b =", params[1])
        dev = np.sqrt(np.diag(cov))
        print("  | a dev =", dev[0])
        print("  | b dev =", dev[1])
        y_fit = lorenzian(x, params[0], params[1])
        print("  | RMSE =", np.linalg.norm(y-y_fit))
        return params
    except RuntimeError:
        return None


# def lorenzian_asym(w, a1, b1, a2, b2):
#     return -1./2. * ( a1 / (1j*w - b1) - a2 / (1j*w + b2) )
def lorenzian_asym(params, x):
    a1 = params[0]
    b1 = params[1]
    a2 = params[2]
    b2 = params[3]
    return -1./2. * ( a1 / (1j*x - b1) - a2 / (1j*x + b2) )


def func_residual(params, x, y):
    return abs( y - lorenzian_asym(params,x) )


def fit_asym(x, z):
    params_sym=fit_sym(x, z.real)
    if params_sym is None:
        return None

    print("  fit with asymmetric function")
    try:
        # params, cov = curve_fit(lorenzian_asym, w, basis)
        init_params = np.array(list(params_sym) * 2)
        params, cov = leastsq(func_residual, init_params, args=(x,z))
        if params[0]<0:
            # fix the sign of parameters so that a1>0, using the transformation
            # a1, b1, a2, b2 --> -a2, -b2, -a1, -b1
            params = [ -params[2], -params[3], -params[0], -params[1] ]
        print("  | a1 =", params[0])
        print("  | b1 =", params[1])
        print("  | a2 =", params[2])
        print("  | b2 =", params[3])
        z_fit = np.array([ lorenzian_asym(params, v) for v in x ])
        print("  | RMSE =", np.linalg.norm(z-z_fit))
        return z_fit
    except RuntimeError:
        return None


def plot_common(filename):
    plt.legend()
    # plt.xlim([-30,30])
    plt.xlabel("$\omega_n$")
    plt.axhline(0) # y=0
    plt.savefig(filename)
    print(" '%s'" %filename)
    plt.close()


def plot_common_3d(xy, z, filename):
    x, y = np.meshgrid(xy, xy)
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.plot_surface(x, y, z, cmap=cm.PuBu)
    plt.xlabel("$\omega_n$")
    plt.ylabel("$\omega_{n'}$")
    plt.savefig(filename, dpi=200)
    print(" '%s'" %filename)
    plt.close()


def save_common_3d(xy, z, filename):
    """Save data for 3D plot using gnuplot"""
    with open(filename, 'w') as f:
        for i,x in enumerate(xy):
            for j,y in enumerate(xy):
                print(x, y, z[i,j], file=f)
            print("", file=f)
    print(" '%s'" %filename)


def print_warning(_str):
    framerecords = inspect.stack()
    # print framerecords[1]
    _, filename = os.path.split( framerecords[1][1] )
    methodname = framerecords[1][3]
    print("*******************************************************")
    print("*** WARNING: in %s() in file '%s'" %(methodname, filename))
    print("***", _str)
    print("*******************************************************")


def sum_w(_u, _beta):
    assert isinstance(_u, np.ndarray)
    assert _u.ndim == 1

    # TODO: 1/iw tail in u(iw)
    u_sum = np.sum(_u) / _beta
    return u_sum


class G2SCL_core:
    """
    Find decoupling of G2
        X(iw,iw') = Phi_L(iw) * Phi_R(iw')

    First call one of set_* functions, and then call get_* functions

    set functions
        set_exact     : SVD of full X(iw,iw')
        set_anti_diag : X(iw,-iw), Phi_L becomes real
        set_1dfull    : X(iw,iw_n) with fixed n
        set_1dhalf    : [recommended] X(iw,-iw_n) for w>0 with fixed n

    get functions
        s0     : [float] The largest singular value
        Phi_L  : [1D ndarray] Phi_L(iw)
        u0     : [1D ndarray] u0(iw)
        singularvalues : [1D ndarray] SV[i]
    """

    def __init__(self, beta=math.pi, handle_exception='zero'):
        """
        beta for plot, 'use_chiloc' option in get_Phi_L(), and set_2pole()

        handle_exception [string]:
            specifies behaviors when exception is raised
                'zero' : replace invalid results with zero
                'ignore' : use invalid results neglecting exceptions
                'exit' : exit immediately
        """
        self.beta = beta
        self.valid = True

        assert handle_exception in ['zero', 'ignore', 'exit']
        self.handle_exception = handle_exception

        # quantities to be set
        self.s0 = None
        self.u0 = None
        self.v0c = None
        self.singular_values = None
        self.X_exact = None
        self.chi_loc = None
        self.u0_svd = None  # TODO: obsolete?

    # =========================================================================
    # set functions

    def set_exact(self, _x):
        """
        Input
            _x     : [2D ndarray] X(iw,iw')
        """
        isinstance(_x, np.ndarray)
        assert _x.ndim == 2

        self.X_exact = _x

        mat_u, vec_s, mat_vh = np.linalg.svd(_x)
        u0 = mat_u[:,0]  # get the first column vector as an array
        v0c = mat_vh[0,:]  # get the first row vector as an array

        # shift phase so that <u0> and <v0> be real
        z = np.sum(u0)
        u0 = shift_phase(u0, z.conj())
        v0c = shift_phase(v0c, z)
        self.u0 = u0
        self.v0c = v0c
        self.s0 = vec_s[0]
        self.singular_values = vec_s.copy()

        # NOTE
        # assert np.allclose(u0, u0[::-1].conj())    # u0(iw) = u0(-iw)^*
        # assert np.allclose(v0c.conj(), v0c[::-1])  # v0(iw) = v0(-iw)^*
        # assert np.allclose(u0, v0c[::-1].conj())   # u0(iw) = v0(-iw)
        try:
            assert np.allclose(u0, u0[::-1].conj()),   "u0(iw) = u0(-iw)^* is not fulfilled"
            assert np.allclose(v0c.conj(), v0c[::-1]), "v0(iw) = v0(-iw)^* is not fulfilled"
            # NOTE: may not be fulfilled in CTQMC
            # TODO
            assert np.allclose(u0, v0c[::-1].conj()),  "u0(iw) = v0(-iw) is not fulfilled"
        except AssertionError as e:
            print_warning(e)
            self.valid = False
            if self.handle_exception == 'exit':
                exit(1)

    def set_anti_diag(self, _x, flag_fit=False):
        """
        Input
            _x     : [1D ndarray] X(iw,-iw) for w>0
        """

        isinstance(_x, np.ndarray)
        assert _x.ndim == 1

        # print _X
        # for x in _x:
        #     # TODO
        #     # assert x.imag < 1e-10  # NOTE: may not be fulfilled in CTQMC
        #     pass

        # NOTE
        try:
            for x in _x:
                assert x.real >= 0, "X >= 0 is not fulfilled"
        except AssertionError as e:
            print_warning(e)
            self.valid = False
            if self.handle_exception == 'exit':
                exit(1)

        if flag_fit:
            pass
        else:
            # phi = np.sqrt(_X.real)  # element-wise square root
            phi = np.sqrt(abs(_x.real))  # element-wise square root
            phi = np.append(phi[::-1].conj(), phi)  # extend phi to negative frequency side
            # s0 = np.linalg.norm(phi)
            # s0 = np.linalg.sum(phi*phi)
            self.s0 = np.dot(phi,phi)
            self.u0 = phi / math.sqrt(self.s0)

        # self.__calc_v0c()
        self.__calc_v0c( +1 if self.valid else -1 )

    def set_1dfull(self, _x, n, flag_fit=False):
        """
        Input
            _x     : [1D ndarray] X(iw,iw_n)
            n      : [int] Matsubara index (0 means w_0 = pi*T)
        """

        isinstance(_x, np.ndarray)
        assert _x.ndim == 1

        w0 = _x.shape[0] / 2
        x_diag = _x[w0 - 1 - n]  # X(-iw_n, iw_n)
        assert x_diag.imag < 1e-10
        assert x_diag.real > 0

        if flag_fit:
            pass
        else:
            phi_n_abs = np.sqrt(x_diag.real)
            phi = _x / phi_n_abs
            # fix phase
            fac = phase( np.sum(phi) )
            phi *= fac.conj()

            # s0 = np.linalg.norm(phi)
            # s0 = np.linalg.sum(phi*phi)
            self.s0 = np.dot(phi.conj(),phi).real
            self.u0 = phi / math.sqrt(self.s0)

        self.__calc_v0c()

    def set_1dhalf(self, _x, n, z, flag_fit=False):
        """
        Input
            _x : [1D ndarray] X(iw,-iw_n) for w>0
            n  : [int] Matsubara index (0 means w_0 = pi*T). n>=0 (w_n>0).
            z  : [complex] X(iw_m,-iw_n) + X(-iw_m,-iw_n)
                  m is arbitrary, but m=n should be avoided.
                  The phase of phi is determined so that z becomes real.
        """

        isinstance(_x, np.ndarray)
        assert _x.ndim == 1
        assert n >= 0

        x_diag = _x[n]
        # TODO
        # assert x_diag.imag < 1e-10  # NOTE: may not be fulfilled in CTQMC

        # NOTE
        # assert x_diag.real > 0
        try:
            assert x_diag.real > 0, "x_diag > 0 is not fulfilled"
        except AssertionError as e:
            print_warning(e)
            self.valid = False
            if self.handle_exception == 'exit':
                exit(1)

        if flag_fit:
            pass
        else:
            phi_n_abs = np.sqrt(abs(x_diag.real))
            phi = _x / phi_n_abs
            # fix phase
            phi *= phase(z).conj()
            phi = np.append(phi[::-1].conj(), phi)  # extend phi to negative frequency side

            # s0 = np.linalg.norm(phi)
            # s0 = np.linalg.sum(phi*phi)
            self.s0 = np.dot(phi.conj(),phi).real
            self.u0 = phi / math.sqrt(self.s0)

        # self.__calc_v0c()
        self.__calc_v0c( np.sign(x_diag.real) )

    def set_2pole(self, energy_plus, energy_minus, nw):
        """
        Input
            E_plus  : [real] Excitation energy to n+1 state (>0)
            E_minus : [real] Excitation energy to n-1 state (>0)
            N_w     : [int] Cutoff of Matsubara frequency
        """

        assert energy_plus > 0
        assert energy_minus > 0
        assert nw % 2 == 0

        # Matsubara frequency
        iw = np.array([(2*n+1) * math.pi / self.beta * 1.j for n in range(-nw // 2, nw // 2)])

        # phi(iw) : normalized by  T*sum_w phi(iw) = 1
        phi = 1. / (iw + energy_minus) - 1. / (iw - energy_plus)

        # u0(iw) : normalized by  sum_w | u0(iw) |^2 = 1
        self.u0 = phi / np.linalg.norm(phi)

        # use_chiloc=True option should be specified to use s0_approx instead of s0
        self.s0 = 0

        self.__calc_v0c()
        # self.__calc_v0c( np.sign(x_diag.real) )

    def __calc_v0c(self, sign=1):
        # u0(iw) = v0(-iw) * sign
        self.v0c = self.u0[::-1].conj() * sign

    def set_chiloc(self, chiloc):
        self.chi_loc = chiloc

    # =========================================================================
    # get functions

    def get_Phi_L(self, use_chiloc=False):
        """
        If use_chiloc==True, s0 is replaced with the approximated value computed from chi_loc.
        This mode is used for computing chi(q) in the RPA-type SCL formula.
        """
        if self.valid or self.handle_exception == 'ignore':
            # return math.sqrt(self.s0) * self.u0
            s0 = self.get_s0(use_chiloc)
            return math.sqrt(s0) * self.u0
        elif self.handle_exception == 'zero':
            return np.zeros(self.u0.shape)
        else:
            raise Exception("should not reach here")

    def get_s0(self, use_chiloc=False):
        if not use_chiloc:
            return self.s0
        else:
            return self.__get_s0_approx()

    def __get_s0_approx(self):
        """
        Return the approximated value of s0 computed from chi_loc
        for checking validity of the approximation
        """

        # calculate s0_approx if not calculated
        if not hasattr(self, 's0_approx'):
            self.s0_approx = self.chi_loc / ( self.beta * sum_w(self.u0, self.beta) * sum_w(self.v0c, self.beta) )

        try:
            assert abs(self.s0_approx.imag) < 1e-10, "s0_approx is not real:" + str(self.s0_approx)
            assert self.s0_approx.real >= 0, "s0_approx is not positive:" + str(self.s0_approx)
        except AssertionError as e:
            print_warning(e)
            self.valid = False
            self.s0_approx = 0
            if self.handle_exception == 'exit':
                exit(1)

        self.s0_approx = self.s0_approx.real
        return self.s0_approx

    def get_u0(self):
        return self.u0

    def get_singularvalues(self):
        return self.singular_values

    # =========================================================================
    # plot functions

    def plot_u0(self, filename, fit=True, data_out=None):
        w=matsubara_mesh(self.u0.real, self.beta)
        plt.plot(w, self.u0.real, marker='.', label='Re u0')
        plt.plot(w, self.u0.imag, marker='.', label='Im u0')

        u0_fit = None
        if fit:
            # fit_sym(w, self.u0.real)
            u0_fit = fit_asym(w, self.u0)
            if u0_fit is not None:
                plt.plot(w, u0_fit.real, '-', label='fit Re', color='r')
                plt.plot(w, u0_fit.imag, '-', label='fit Im', color='g')

        if self.u0_svd is not None:
            plt.plot(w, self.u0_svd.real, '-o', label='Re u0_{SVD}')
            plt.plot(w, self.u0_svd.imag, '-o', label='Im u0_{SVD}')
        plot_common(filename)

        if isinstance(data_out, str):
            array = [w, self.u0.real, self.u0.imag]
            if u0_fit is not None:
                array.append(u0_fit.real)
                array.append(u0_fit.imag)
            np.savetxt(data_out, np.array(array).T)

    def plot_x(self, _filename, use_chiloc=False, flag_data_out=False):
        w=matsubara_mesh(self.u0.real, self.beta)
        basename, ext = os.path.splitext(_filename)

        x_approx = np.einsum("i,j->ij", self.u0, self.v0c) * self.get_s0(use_chiloc)

        filename = basename + '_approx' + ext
        plot_common_3d(w, x_approx.real, filename)
        if flag_data_out:
            save_common_3d(w, x_approx.real, basename + '_approx.dat')

        if self.X_exact is not None:
            filename = basename + '_exact' + ext
            plot_common_3d(w, self.X_exact.real, filename)
            if flag_data_out:
                save_common_3d(w, self.X_exact.real, basename + '_exact.dat')

            filename = basename + '_diff' + ext
            plot_common_3d(w, self.X_exact.real-x_approx.real, filename)
            if flag_data_out:
                save_common_3d(w, self.X_exact.real-x_approx.real, basename + '_diff.dat')

    def plot_sv(self, _filename, xmax=30):
        if self.singular_values is None:
            return

        basename, ext = os.path.splitext(_filename)
        # print basename, ext

        plt.plot(self.singular_values, '-o')
        plt.legend()
        plt.xlim([0,xmax])
        plt.xlabel("$i$")
        plt.ylabel("$s_i$")
        plt.axhline(0) # y=0
        filename = basename + '-linear' + ext
        plt.savefig(filename)
        print(" '%s'" %filename)

        # semi-log scale
        plt.yscale('log')
        filename = basename + '-log' + ext
        plt.savefig(filename)
        print(" '%s'" %filename)
        plt.close()

    # =========================================================================
    # save functions

    def save_sv(self, filename):
        if self.singular_values is None:
            return
        np.savetxt(filename, self.singular_values)
        print(" '%s'" %filename)
