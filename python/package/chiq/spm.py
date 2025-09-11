import numpy as np
from scipy.integrate import quad
import sys

try:
    import sparse_ir
    SPARSE_IR_AVAILABLE = True
except ImportError:
    SPARSE_IR_AVAILABLE = False

try:
    import admmsolver
    import admmsolver.optimizer
    import admmsolver.objectivefunc
    from admmsolver.objectivefunc import L1Regularizer, ConstrainedLeastSquares, NonNegativePenalty
    from admmsolver.matrix import identity
    from admmsolver.optimizer import SimpleOptimizer

    ADMMSOLVER_AVAILABLE = True
except ImportError:
    ADMMSOLVER_AVAILABLE = False

if SPARSE_IR_AVAILABLE and ADMMSOLVER_AVAILABLE:
    SPM_AVAILABLE = True
else:
    SPM_AVAILABLE = False

class SpM:
    def __init__(self, beta, wmax, matsubara_points, chi_iwn):
        can_use = True
        if not SPARSE_IR_AVAILABLE:
            print("ERROR: sparse_ir is not installed.")
            can_use = False
        if not ADMMSOLVER_AVAILABLE:
            print("ERROR: admmsolver is not installed.")
            can_use = False
        if not can_use:
            print("ERROR: SpM method is not available.")
            sys.exit(1)

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
        self.A = np.einsum("il,l->il", U, S)
        y = self.g_tau

        C = (self.A[0] + self.A[-1]).reshape(1, -1)
        D = np.array([y[0] + y[-1]])

        self.lstsq_F = ConstrainedLeastSquares(0.5, A=self.A, y=y, C=C, D=D)


    def fit(self, l1_coeff, max_iter=1000, initial_mu=1.0):
        nw = 101
        l1_F = L1Regularizer(l1_coeff, self.basis.size)
        nonneg_F = NonNegativePenalty(nw)
        ws = np.linspace(0.0, self.wmax, nw)
        V = self.basis.v(ws).T

        objective_functions = [self.lstsq_F, l1_F, nonneg_F]
        equality_conditions = [(0, 1, identity(self.basis.size), identity(self.basis.size)),
                               (0, 2, V, identity(nw)),
                               ]
        p = admmsolver.optimizer.Problem(objective_functions, equality_conditions)

        self.opt = SimpleOptimizer(p, mu=initial_mu)
        self.opt.solve(max_iter)


    def training_loss(self):
        return np.sum((self.A @ self.opt.x[0] - self.g_tau) ** 2)


    def optimize_lambda(self, loglambda_min, loglambda_max, loglambda_num, max_iter=1000, initial_mu=1.0):
        loglambdas = np.linspace(loglambda_min, loglambda_max, loglambda_num)
        losses = np.zeros(loglambda_num)
        for i, loglambda in enumerate(loglambdas):
            self.fit(10**loglambda, max_iter=1000, initial_mu=1.0)
            losses[i] = np.log10(self.training_loss())
        scores = (loglambdas - loglambda_min) / (loglambda_max - loglambda_min) * (losses[-1] - losses[0]) + losses[0]
        scores -= losses
        optimal_loglambda = loglambdas[np.argmax(scores)]
        return optimal_loglambda


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
