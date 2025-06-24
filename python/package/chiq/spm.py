import numpy as np
from scipy.integrate import quad
import sparse_ir
import admmsolver
import admmsolver.optimizer
import admmsolver.objectivefunc
from admmsolver.objectivefunc import L1Regularizer, ConstrainedLeastSquares
from admmsolver.matrix import identity
from admmsolver.optimizer import SimpleOptimizer

class SpM:
    def __init__(self, beta, wmax, matsubara_points, chi_iwn, l1_coeff, max_iter=1000, initial_mu=1.0):
        self.beta = beta
        self.wmax = wmax
        self.basis = sparse_ir.FiniteTempBasis("B", beta, wmax)

        for n in matsubara_points:
            assert n % 2 == 0, "Matsubara points must be even"

        sampler_iwn = sparse_ir.MatsubaraSampling(self.basis, matsubara_points, positive_only=True)
        self.g_l = sampler_iwn.fit(chi_iwn)

        sampling_tau = sparse_ir.TauSampling(self.basis)
        self.g_tau = sampling_tau.evaluate(self.g_l)

        ts = sampling_tau.sampling_points
        U = self.basis.u(ts)
        S = self.basis.s
        A = np.einsum("il,l->il", U, S)
        y = self.g_tau

        C = (A[0] + A[-1]).reshape(1, -1)
        D = np.array([y[0] + y[-1]])

        lstsq_F = ConstrainedLeastSquares(0.5, A=A, y=y, C=C, D=D)
        l1_F = L1Regularizer(l1_coeff, self.basis.size)

        objective_functions = [lstsq_F, l1_F]
        equality_conditions = [(0, 1, identity(self.basis.size), identity(self.basis.size))]
        p = admmsolver.optimizer.Problem(objective_functions, equality_conditions)

        self.opt = SimpleOptimizer(p, mu=initial_mu)
        self.opt.solve(max_iter)


    def evaluate(self, ws):

        def A(w):
            V = self.basis.v(w).T
            rho = V @ self.opt.x[0]
            return np.real(rho * np.tanh(0.5*self.beta*w))

        chi_imag = -np.pi * A(ws)
        chi_real = np.zeros_like(ws)
        for iw, w in enumerate(ws):
            chi_real[iw] = quad(A, -self.wmax, self.wmax, weight="cauchy", wvar=w)[0]
        return chi_real + 1j * chi_imag
