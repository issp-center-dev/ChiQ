

import numpy as np
import sys
# import math
# import cmath
from itertools import product, combinations_with_replacement
import copy
import warnings
from collections import namedtuple, OrderedDict, Counter
from fractions import Fraction
from . import point_group_data

decimals = 14
decimals_degene = 11
max_rank = 10
max_num_multipoles_in_one_irrep = 100


def commutator(mat1, mat2):
    return mat1 @ mat2 - mat2 @ mat1

def anticommutator(mat1, mat2):
    return mat1 @ mat2 + mat2 @ mat1

def is_commute(mat1, mat2):
    return np.allclose(mat1 @ mat2, mat2 @ mat1)

def is_anticommute(mat1, mat2):
    return np.allclose(mat1 @ mat2, -mat2 @ mat1)

def is_hermitian(mat):
    # return np.allclose(mat, np.conjugate(mat.T))
    return np.allclose(mat, mat.conj().T)

def is_unitary(mat):
    return np.allclose(mat @ mat.conj().T, np.identity(mat.shape[0]))


def hermitianize(mat):
    return (mat + mat.conj().T) / 2.


# def is_not_duplicate(array):
#     return len(array) == len(set(array))


def make_matrix_zeros(dim):
    return np.zeros((dim, dim), dtype=complex)

def make_matrix_identity(dim):
    return np.identity(dim)

def make_matrix_random_real(dim):
    return np.random.rand(dim, dim)

def make_matrix_random_complex(dim):
    return make_matrix_random_real(dim) + 1.j * make_matrix_random_real(dim)


def isreal(x):
    return np.allclose(x.imag, np.zeros(x.shape))

def matrix_weight(mat):
    # w = np.trace(mat * mat.getH())
    w = matrix_overlap(mat, mat)
    assert isreal(w), f"{w} is not real"
    return w.real

def matrix_normalize(mat):
    w = matrix_weight(mat)
    return mat / np.sqrt(w)

def matrix_overlap(mat1, mat2):
    # return np.trace(mat1 * mat2.getH())
    return np.trace(mat1.conj().T @ mat2)

def matrix_orthogonalize(mat1, mat2):
    # overlap
    # p = matrix_overlap(mat1, mat2)
    p = matrix_overlap(mat2, mat1)
    # print(p)
    r = mat1 - p * mat2
    assert np.isclose(matrix_overlap(r, mat2), 0)
    return r

# def matrix_project(mat1, mat2):
#     matrix_orthogonalize(mat1, mat2)


def make_rot_mat(generator, theta):
    """
    compute R = e^{i theta generator}

    input:
        generator  np.ndarray
        theta  float
    return:
        R  np.ndarray
    """
    assert isinstance(generator, np.ndarray)

    # shape = generator.shape
    # return np.asmatrix(np.zeros(shape))

    W, V = np.linalg.eigh(generator)
    # V = np.asmatrix(V)

    # check V^h * G * V = W
    # print("V^h.G.V = ", V.getH() * generator * V)
    assert np.allclose(V.conj().T @ generator @ V, np.diag(W))

    # V * e^{iWt} * V^h
    # e_iwt = [cmath.exp(1.j * theta * w) for w in W]
    # e_W = np.asmatrix(np.diag(e_iwt))
    # return V * e_W * V.getH()
    e_iwt = np.exp(1.j * theta * W)
    e_W = np.diag(e_iwt)
    return V @ e_W @ V.conj().T


class Multipole(object):
    def __init__(self, mat, irp=None, rank=None, space=None):
        # self.name = name
        self.irp = irp
        self.rank = rank
        self.space = space
        self.mat = mat

    def __repr__(self):
        return " irp = {}, rank = {}, space = {}, mat.shape = {}".format(
            self.irp, self.rank, self.space, self.mat.shape)


def append_multipoles(multipoles1, multipoles2, func=None, **info):
    """
    apppend multipoles2 to multipoles1
    """
    assert isinstance(multipoles1, OrderedDict)
    assert isinstance(multipoles2, OrderedDict)

    def do_nothing(x):
        return x
    if func is None:
        func = do_nothing

    for (irp1, m1_irp), (irp2, m2_irp) in zip(list(multipoles1.items()), list(multipoles2.items())):
        assert irp1 == irp2  # ensured by OrderedDict
        for m in m2_irp:
            m_info = m.__dict__
            mat = m_info.pop('mat')
            m_info.update(info)
            m1_irp.append(Multipole(mat=func(mat), **m_info))


class PointGroup(object):

    # Operator = namedtuple('Operator', ('name', 'mat'))
    EigenStates = namedtuple('EigenStates', ('name', 'dim', 'vecs'))

    def __init__(self, angular_momentum, group, basis_order=None, verbose=True):

        self.verbose = verbose

        # 2J as int
        twoj = Fraction(angular_momentum) * 2
        assert twoj.denominator == 1, "J must be integer or half-integer"
        self.twoj = twoj.numerator
        self.flag_half_integer = self.twoj % 2 == 1

        # self.dim = 2*l+1
        self.dim = self.twoj+1

        # self._make_angular_momentum_matrices(l)
        self._make_angular_momentum_matrices(self.twoj, basis_order)
        self._set_point_group(group)
        self._make_symmetry_operation_matrices()
        self.restore_states()  # set self._reduce_states_vecs, self._reduce_states_space

    @staticmethod
    def group_list():
        return list(point_group_data.PointGroupData.keys())

    def _print_verbose(self, *args):
        if self.verbose:
            print(*args)

    def _make_angular_momentum_matrices(self, twoj, basis_order):
        # operator Lz
        # self.Lz = np.matrix(np.diag([lz for lz in range(-l, l + 1)]))
        # self.Lz = np.matrix(np.diag([twojz/2. for twojz in range(-twoj, twoj + 2, 2)]))
        self.Lz = np.diag([twojz/2. for twojz in range(-twoj, twoj + 2, 2)])

        # operator L+
        Lp = make_matrix_zeros(self.dim)
        # for lz in range(-l, l):
        #     i = lz+l
        for i, twojz in enumerate(range(-twoj, twoj, 2)):
            lz = twojz / 2.
            l = twoj / 2.
            Lp[i + 1, i] = np.sqrt((l - lz) * (l + lz + 1))

        # operator L-
        Lm = Lp.T

        self.Lx = (Lp + Lm) / 2
        self.Ly = (Lp - Lm) / 2j

        # check [Lz,L+] = L+
        assert np.allclose(commutator(self.Lz, Lp), Lp)
        # check [Lz,L-] = L-
        assert np.allclose(commutator(self.Lz, Lm), -Lm)
        # check [L+,L-] = 2 Lz
        assert np.allclose(commutator(Lp, Lm), 2 * self.Lz)

        if basis_order is not None:
            assert isinstance(basis_order, (list, tuple, np.ndarray))
            basis_order = np.asarray(basis_order)
            assert basis_order.shape == (self.dim,), f"basis_order={basis_order}, dim={self.dim}"
            self.Lx = self.Lx[np.ix_(basis_order, basis_order)]
            self.Ly = self.Ly[np.ix_(basis_order, basis_order)]
            self.Lz = self.Lz[np.ix_(basis_order, basis_order)]

    @property
    def irreps(self):
        return [irp.name for irp in self.G.irreps1 + self.G.irreps2]

    @property
    def irreps1(self):
        return [irp.name for irp in self.G.irreps1]

    @property
    def basis(self):
        return np.diag(self.Lz)

    def _set_point_group(self, group):
        self._print_verbose("\nimporting data of point group %s..." % group)
        self.G = point_group_data.PointGroupData[group](double_group=self.flag_half_integer)

    def _make_symmetry_operation_matrices(self):
        self._print_verbose("\ncreating matrices of symmetry operations...")

        for op in self.G.operators1 + self.G.operators2:
            self._print_verbose(op.name, op.optype, op.axis, op.angle)

            if op.optype == 'identity':
                mat = make_matrix_identity(self.dim)
            elif op.optype == 'rot':
                assert isinstance(op.angle, str)
                angle = Fraction(op.angle)
                theta = 2 * np.pi * float(angle)

                # print(op.axis)
                assert isinstance(op.axis, str)
                xyz = op.axis.split()
                # print(xyz)

                # for syntax like "x, y, -z"
                if all(x in ['x', 'y', 'z', '-x', '-y', '-z'] for x in xyz):
                    assert len(xyz) <= 3

                    gen_ops = {'x': self.Lx, 'y': self.Ly, 'z': self.Lz,
                               '-x': -self.Lx, '-y': -self.Ly, '-z': -self.Lz}
                    gen = make_matrix_zeros(self.dim)
                    for x in xyz:
                        gen += gen_ops[x]
                    gen /= np.sqrt(len(xyz))

                # for syntax like "1/6 xy"
                elif xyz[1] in ['xy', 'yz', 'zx'] and len(xyz) == 2:
                    assert len(xyz) == 2

                    plane = xyz[1]
                    phi = 2 * np.pi * Fraction(xyz[0])

                    gen_ops = {'xy': [self.Lx, self.Ly], 'yz': [self.Ly, self.Lz], 'zx': [self.Lz, self.Lx]}
                    gen = np.cos(phi) * gen_ops[plane][0] + np.sin(phi) * gen_ops[plane][1]
                else:
                    raise ValueError

                mat = make_rot_mat(gen, theta)
            elif op.optype == 'mirror':
                raise NotImplementedError
            else:
                raise ValueError('op.optype =', op.optype)

            if op.inversion:
                if not self.flag_half_integer:
                    # if self.l % 2:
                    if self.twoj % 4:
                        mat *= -1.
                else:
                    raise NotImplementedError
                    # TODO: for J with half-integer

            if op.double:
                mat *= -1.

            op.mat = mat

        # Time reversal operator
        self.T_mat = make_rot_mat(self.Ly, np.pi)

        # store data for reduce_matrix
        self.dim_full = self.dim
        # self.Ops_full = copy.deepcopy(self.Ops)
        for op in self.G.operators:
            op.backup()

    def weights_char(self, chars):
        # weights = self.projection(chars)
        weights = []
        for _chars in self.G.character_table:
            w = self.G.inner_product(chars, _chars) / self.G.group_order
            # print(p)
            # weights.append(np.round(w, decimals))
            weights.append(w)
        # for irp, w in zip(self.G.irreps, weights):
        #     print(irp.name, w)
        return weights

    def _decompose_mat(self, mat, time_reversal):
        # assert ishermitian(mat)

        # print("original matrix:")
        # print(mat)

        mat_projected = []
        # for irp in Irps:
        # for irp, chars in zip(self.G.irreps, self.G.character_table):
        for irp in self.G.irreps:
            # print(irp.name, irp.dim)
            # print(chars)

            # mat_p = np.asmatrix(np.zeros((self.dim, self.dim), dtype=complex))
            mat_p = make_matrix_zeros(self.dim)
            for op in self.G.operators:
                mat_p += op.mat @ mat @ op.mat.conj().T * self.G.characters_dict[(irp.name, op.name)]
            # mat_p /= self.G.group_order
            mat_p *= irp.dim / float(self.G.group_order)
            mat_p = np.round(mat_p, decimals)
            # print(mat_p)
            # assert is_hermitian(mat_p)

            # Time reversal symmetry
            if time_reversal is not None:
                if time_reversal == 'even':
                    s = +1
                elif time_reversal == 'odd':
                    s = -1
                else:
                    raise ValueError(f"Unexpected value: time_reversal = {time_reversal}")
                mat_p += s * self.T_mat @ np.conj(mat_p) @ self.T_mat.conj().T
                mat_p /= 2

            mat_projected.append(mat_p)
        return mat_projected

    def decompose_mat(self, mat, time_reversal=False):
        assert isinstance(mat, np.ndarray)
        assert mat.shape == (self.dim, self.dim)

        if time_reversal:
            mat_projected = self._decompose_mat(mat, time_reversal='even')
            mat_projected += self._decompose_mat(mat, time_reversal='odd')
        else:
            mat_projected = self._decompose_mat(mat, time_reversal=None)

        # check if weight is not lost
        mat_reconstructed = make_matrix_zeros(self.dim)
        for mat_p in mat_projected:
            mat_reconstructed += mat_p
        assert np.allclose(mat, mat_reconstructed)

        return mat_projected

    def weights_mat(self, mat, time_reversal=False):
        """Decompose matrix into irreducible representations and return those weights

        Args:
            mat (np.ndarray): matrix

        Returns:
            list(float): weights[i_irp]
        """
        assert isinstance(mat, np.ndarray)
        assert mat.shape == (self.dim, self.dim)

        mat_projected = self.decompose_mat(mat, time_reversal)

        weights = []
        for mat_p in mat_projected:
            weights.append(matrix_weight(mat_p))

        if time_reversal:
            assert len(weights) == len(self.G.irreps) * 2  # for even, odd
        else:
            assert len(weights) == len(self.G.irreps)

        return weights

    def weights_mat_multipoles(self, mat):
        """Decompose matrix into irreducible representations and rank and return those weights

        Args:
            mat (np.ndarray): matrix

        Returns:
            OrderedDict(list(float)): weights[irp][rank]
        """
        assert isinstance(mat, np.ndarray)
        assert mat.shape == (self.dim, self.dim)

        multipoles = self.multipoles

        # weights = OrderedDict()
        weights = self._multipoles_empty()
        w_sum = 0
        mat_recovered = make_matrix_zeros(self.dim)
        for irp, m_irp in list(multipoles.items()):
            weights[irp] = [0] * max_rank
            # print("***", irp)
            for m in m_irp:
                # print(matrix_weight(m.mat))
                coef = matrix_overlap(mat, m.mat) / matrix_weight(m.mat)
                w = (coef * coef.conj()).real
                weights[irp][m.rank] += w
                # print(np.linalg.norm(w))
                # print(w)
                w_sum += w
                mat_recovered += coef * m.mat
        assert np.isclose(w_sum, matrix_weight(mat))
        # assert np.allclose(mat, mat_recovered)
        if not np.allclose(mat, mat_recovered, atol=1e-3):
            print("Warning: Decomposition into rank is not prefect.", file=sys.stderr)
        return weights

    def _gen_rank_operators(self, rank):
        """
        generator of operators in a given rank
        """
        if rank == 0:
            yield make_matrix_identity(self.dim_full)
        else:
            for op_low in self._gen_rank_operators(rank-1):
                for op in [self.Lx, self.Ly, self.Lz]:
                    yield op_low @ op
            # TODO: store results

    # def _gen_rank_operators_core(self, rank):
    #     """
    #     generator of operators in a given rank
    #     """
    #     import functools
    #     import operator
    #
    #     if rank == 0:
    #         yield make_matrix_identity(self.dim_full)
    #     else:
    #         # for comb in combinations_with_replacement([self.Lx, self.Ly, self.Lz], rank):
    #         #     yield functools.reduce(operator.mul, comb)
    #         # for ops in product([self.Lx, self.Ly, self.Lz], repeat=rank):
    #         #     yield functools.reduce(operator.mul, ops)
    #         for op_low in self._gen_rank_operators(rank-1):
    #             for op in [self.Lx, self.Ly, self.Lz]:
    #                 yield op_low * op
    #         # TODO: store results
    #
    # def _gen_rank_operators(self, rank):
    #     if not hasattr(self, '_store_rank_operators'):
    #         self._store_rank_operator = [None] * max_rank
    #
    #     if self._store_rank_operator[rank] is None:
    #         self._store_rank_operator[rank] = [op for op in self._gen_rank_operators_core(rank)]
    #
    #     for op in self._store_rank_operator[rank]:
    #         yield op

    def calc_chars(self, vecs, integer=False):
        assert vecs.shape[0] == self.dim

        chars = []
        syms = []
        for op in self.G.operators:
            char = np.trace(vecs.conj().T @ op.mat @ vecs)
            char = np.round(char, decimals_degene)
            # print(op.name, char)
            if op.name not in syms:
                chars.append(char)
                syms.append(op.name)
            else:
                assert np.isclose(chars[-1], char)

        if integer:
            # TODO
            pass

        return chars

    @property
    def cf_states(self):
        """
        return:
            list of namedtuple('EigenStates', ('name', 'dim', 'vecs'))
        """
        if not hasattr(self, '_cf_states'):
            # if self.flag_half_integer:
            #     self.G.switch_double()  # switch to double rep

            self._cf_states = self._make_crystal_field_eigenstates()
            self._recombine_cf_states()
            # TODO: check eigenstates
            # TODO: fix the phase

            for cf in self._cf_states:
                assert cf.vecs.shape == (self.dim, cf.dim)  # column vectors
        return self._cf_states

    @property
    def unitary_matrix(self):
        """
        return:
            np.ndarray  < lz | cf >
        """
        cf_states = self.cf_states
        mat = np.hstack([cf.vecs for cf in cf_states])
        assert mat.shape == (self.dim, self.dim)
        return mat

    def convert_to_cf_states(self, mat):
        """
        return:
            np.ndarray  < irrep | mat | irrep' >
        """
        isinstance(mat, np.ndarray)

        unimat = self.unitary_matrix
        # return unimat * mat * unimat.getH()
        return unimat.conj().T @ mat @ unimat

    def _make_crystal_field_eigenstates(self):
        randmat_r = make_matrix_random_real(self.dim)
        randmat_r = hermitianize(randmat_r)

        # point-charge potential <- random matrix within A1g irrep
        # TODO: simplify
        self.G.use_single_reps()
        pot = self.decompose_mat(randmat_r)[0]  # A1g
        if self.flag_half_integer:
            self.G.use_double_reps()  # switch to double rep

        eigval, eigvec = np.linalg.eigh(pot)
        # eigvec = np.round(eigvec, decimals)
        # print(" eigenvalues =", eigval)
        # print(eigvec)

        eigval = np.round(eigval, decimals_degene)
        assert is_unitary(eigvec)

        # check degeneracy
        c = Counter(eigval)
        self._print_verbose(c)

        # classify eigenvectors
        dims = []
        eigvec_split = []
        offset = 0
        for x in eigval:
            # deg = c[x]
            deg = c.pop(x, None)
            if deg is not None:
                eigvec_split.append(eigvec[:, offset:offset+deg])
                dims.append(deg)
                offset += deg
                # print(deg)
        self._print_verbose(" dims = {}".format(dims))
        # print("num irrep = ", len(eigvec_split))

        eig_irreps = []

        for vecs in eigvec_split:
            # print(vecs.shape)

            # compute characters
            chars = self.calc_chars(vecs, integer=True)
            # print(chars)

            # weights of each irrep
            weights = self.weights_char(chars)
            # print(weights)

            # check if irrep is uniquely identified:
            #   weights should be something like [0, ..., 0, 1, 0, ..., 0]
            def one_or_zero(_w):
                nums = (0, 1)
                for n in nums:
                    if np.isclose(_w, n):
                        return n
                raise Exception
            weights_10 = [one_or_zero(w) for w in weights]

            # check if 1 appears only once
            # import pdb; pdb.set_trace()
            assert sum(weights_10) == 1

            pos = weights_10.index(1)
            eig_irreps.append(self.G.irreps[pos].name)

        self._print_verbose(" irreps = {}".format(eig_irreps))

        # store data into namedtuple
        assert len(eig_irreps) == len(dims) == len(eigvec_split)
        eigens = [PointGroup.EigenStates(eig_irreps[i], dims[i], eigvec_split[i]) for i in range(len(dims))]
        # sort by irrep name
        eigens.sort(key=lambda eigen: self.irreps.index(eigen.name))
        return eigens

    def reduce_states(self, cf_list):
        """
        input:
            vecs  np.matrix or np.ndarray
                  with shape (dim, <=dim)
            cf_list  PointGroup.EigenStates, or list or tuple of PointGroup.EigenStates
        """
        if isinstance(cf_list, PointGroup.EigenStates):
            vecs = cf_list.vecs.copy()
            space = (cf_list.name, )
        else:
            vecs = np.hstack([cf.vecs for cf in cf_list])
            space = tuple([cf.name for cf in cf_list])

        assert vecs.shape[0] == self.dim_full
        assert vecs.shape[1] <= self.dim_full

        # TODO: check orthogonality
        # TODO: check normalization

        self.dim = vecs.shape[1]
        self._reduce_states_vecs = vecs
        self._reduce_states_space = space

        def reduce_matrix(mat):
            return vecs.conj().T @ mat @ vecs

        # make operators defined in the reduced space
        for op in self.G.operators:
            op.backup()
            op.mat = reduce_matrix(op.mat)

        self.mat_cf_reduce = vecs.copy()

        # reset
        if hasattr(self, '_multipoles'):
            del self._multipoles
        if hasattr(self, '_multipoles_labeled'):
            del self._multipoles_labeled
        if hasattr(self, '_multipoles_original_basis'):
            del self._multipoles_original_basis

    def restore_states(self):
        # restore operators
        for op in self.G.operators:
            op.restore()
        # restore states
        self.reduce_states(PointGroup.EigenStates('full', self.dim_full, make_matrix_identity(self.dim_full)))
        self.mat_cf_reduce = np.identity(self.dim_full)

    # def project_mat_onto_cf_states(self, mat, vecs):
    #     assert isinstance(mat, np.matrix)
    #     assert isinstance(vecs, np.matrix)
    #     assert mat.shape[0] == mat.shape[1] == self.dim_full
    #     assert vecs.shape[0] == self.dim_full
    #
    #     assert mat.shape
    #     proj = vecs * vecs.getH()
    #     return proj * mat * proj

    @property
    def multipoles(self):
        """
        return:
            dict of list of class Multipole
        """
        if not hasattr(self, '_multipoles'):
            self.G.use_single_reps()  # switch to single rep
            # self._multipoles = self._make_multipole_operators()
            self._multipoles = self._make_multipole_operators_ranked()
            assert np.sum([len(m_irp) for m_irp in list(self._multipoles.values())]) == self.dim ** 2
            self._recombine_multipoles()
            self._check_multipole_orthogonality(self._multipoles)
            # TODO: how to check completeness?
        return self._multipoles

    @property
    def multipoles_original_space(self):
        """
        return:
            dict of list of class Multipole
        """
        if not hasattr(self, '_multipoles_original_basis'):
            self._multipoles_original_basis = copy.deepcopy(self._multipoles)
            for m_irp in list(self._multipoles_original_basis.values()):
                for m in m_irp:
                    m.mat = self.mat_cf_reduce @ m.mat @ self.mat_cf_reduce.conj().T
        return self._multipoles_original_basis

    @property
    def multipoles_labeled(self):
        """
        return:
            dict of list of multipole operators

        NOTE: may not be necessary
        """
        if not hasattr(self, '_multipoles_labeled'):
            self._multipoles_labeled = self._make_multipole_operators_labeled()
            assert len(self._multipoles_labeled) == self.dim_full ** 2
            self._check_multipole_orthogonality(self._multipoles_labeled)
        return self._multipoles_labeled

    def get_multipoles(self, irp=None, rank=None, space=None):
        m_hit = self.multipoles_labeled
        if irp is not None:
            m_hit = [m for m in m_hit if m.irp == irp]
        if rank is not None:
            m_hit = [m for m in m_hit if m.rank == rank]
        if space is not None:
            m_hit = [m for m in m_hit if m.space == space]
        return m_hit

    def _multipoles_empty(self):
        multipoles = OrderedDict()
        for irp in self.G.irreps:
            multipoles[irp.name] = []
        return multipoles

    def _check_multipole_orthogonality(self, multipole_list_or_dict):
        self._print_verbose("\nchecking orthogonality of multipole operators...")

        if isinstance(multipole_list_or_dict, list):
            m_list = multipole_list_or_dict
        elif isinstance(multipole_list_or_dict, dict):
            m_list = [m for m_irp in list(multipole_list_or_dict.values()) for m in m_irp]
        else:
            raise Exception("input must be list or dict, but {}".format(type(multipole_list_or_dict)))

        for (i, m1), (j, m2) in product(enumerate(m_list), repeat=2):
            # print(i, m1)
            # print(j, m2)
            # ol = np.trace(m1.mat * m2.mat.getH())
            ol = matrix_overlap(m1.mat, m2.mat)
            # print(ol)
            if i == j:
                assert np.isclose(ol, 1), "ol -1 = {}".format(ol-1)
            else:
                assert np.isclose(ol, 0), "ol = {}".format(ol)

    def _make_multipole_operators(self, generator=None, exclude=None):
        self._print_verbose("_make_multipole_operators")

        multipoles = self._multipoles_empty()
        assert isinstance(multipoles, dict)

        if exclude is None:
            exclude = self._multipoles_empty()

        for n in range(max_num_multipoles_in_one_irrep):
            # randmat = make_matrix_random_complex(self.dim)
            if generator is None:
                randmat = make_matrix_random_complex(self.dim)
            else:
                randmat = generator()
            mat_decomposed = self.decompose_mat(randmat)

            flag_finish = True
            for mat, irp in zip(mat_decomposed, self.G.irreps):
                # print(irp.name)
                # orthogonalize with already obtained multipoles
                for op_obtained in multipoles[irp.name] + exclude[irp.name]:
                    mat = matrix_orthogonalize(mat, op_obtained.mat)
                # append mat if its weight is non-zero
                if not np.isclose(matrix_weight(mat), 0):
                    multipoles[irp.name].append(
                        Multipole(mat=matrix_normalize(mat), irp=irp.name, space=self._reduce_states_space)
                    )
                    flag_finish = False
            if flag_finish:
                break

        for irp in self.G.irreps:
            self._print_verbose(" {}  {}".format(irp.name, len(multipoles[irp.name])))

        return multipoles

    def _make_randmat_rank(self, rank):
        """make random matrix within a given rank"""

        # return a random number in the range [-1:1)
        def rand_abs1():
            x = np.random.rand()  # [0:1)
            return 2. * x - 1

        mat = make_matrix_zeros(self.dim_full)
        for mat_rank in self._gen_rank_operators(rank):
            mat += mat_rank * rand_abs1()

        # return mat
        return self._reduce_states_vecs.conj().T @ mat @ self._reduce_states_vecs

    def _make_multipole_operators_ranked(self):
        multipoles = self._multipoles_empty()
        assert isinstance(multipoles, dict)
        count = 0
        for rank in range(max_rank):
            # print("rank =", rank)
            multipoles_new = self._make_multipole_operators(
                generator=lambda: self._make_randmat_rank(rank),
                exclude=multipoles,
            )
            count += np.sum([len(m_irp) for m_irp in list(multipoles_new.values())])
            append_multipoles(multipoles, multipoles_new, func=None, rank=rank)
            if count >= self.dim**2:
                break
        return multipoles

    def _make_multipole_operators_labeled(self):
        """
        NOTE: may not be necessary
        """

        # make random matrix within a given rank
        def randmat_oddfiagonal(_rank, _offdiagonal=None):
            mat = self._make_randmat_rank(_rank)
            # mat = _vecs.getH() * mat * _vecs

            # if offdiagonal is given, leave only offdiagonal block
            if _offdiagonal is not None:
                n1, n2 = _offdiagonal
                # assert n1 + n2 == _vecs.shape[1]
                # delete diagonal block
                mat[0:n1, 0:n1] = 0
                mat[n1:n1 + n2, n1:n1 + n2] = 0

            return mat

        multipoles = self._multipoles_empty()

        for (i, cf1), (j, cf2) in product(enumerate(self.cf_states), repeat=2):
            # print(cf.name)
            if i == j:
                cfs = (cf1,)
                vecs = cf1.vecs
                offdiagonal = None
                # space = (cf1.name,)
            elif i < j:
                cfs = (cf1, cf2)
                vecs = np.hstack([cf1.vecs, cf2.vecs])
                offdiagonal = (cf1.vecs.shape[1], cf2.vecs.shape[1])
                # space = (cf1.name, cf2.name)
            else:
                continue
            self.reduce_states(cfs)

            multipoles_cf = self._multipoles_empty()
            for rank in range(max_rank):
                # print("rank =", rank)
                multipoles_new = self._make_multipole_operators(
                    generator=lambda: randmat_oddfiagonal(rank, offdiagonal),
                    exclude=multipoles_cf,
                )
                append_multipoles(multipoles_cf, multipoles_new, func=None, rank=rank)

            def matrix_expand(mat):
                return vecs @ mat @ vecs.conj().T
            # append_multipoles(multipoles, multipoles_cf, func=matrix_expand, space=space)
            append_multipoles(multipoles, multipoles_cf, func=matrix_expand)

            self.restore_states()

        # flatten dict of list
        return [m for m_irp in list(multipoles.values()) for m in m_irp]

    def _recombine_bases_core(self, bases, func_matrix_element, func_overlap):
        # TODO
        operators_diagonalized_if_degenerate = ['C4^2', 'C4', 'C2x']
        ops = []
        for sym in operators_diagonalized_if_degenerate:
            ops.extend([op.mat for op in self.G.operators if op.name == sym])

        bases_new = []  # list of new column vector
        for op in ops:
            # matrix within a given cf states
            mat = np.zeros((len(bases), len(bases)), dtype=complex)
            for i, j in product(list(range(len(bases))), repeat=2):
                # mat[i, j] = vec[i].getH() * op * vec[j]
                mat[i, j] = func_matrix_element(bases[i], bases[j], op)
            # print(mat)

            # diagonalize
            eigval, eigvec = np.linalg.eig(mat)
            eigval = np.round(eigval, decimals_degene)
            # print(eigval)
            # print(eigvec)

            # loops for unique eigenvalues
            c = Counter(eigval)
            self._print_verbose("eigenvalues of sym op:", c)
            for key, count in list(c.items()):
                if count == 1:
                    col = list(eigval).index(key)
                    # print(col)
                    # v = cf.vecs * eigvec[:, col]
                    # sum_n vec[n] * coef[n]
                    new_basis = sum([bases[n] * eigvec[n, col] for n in range(len(bases))])
                    new_basis = np.round(new_basis, decimals)

                    # check if the vector v is new
                    flag_new = True
                    for b in bases_new:
                        # print("overlap =", func_overlap(new_basis, b))
                        # if not np.isclose(matrix_overlap(v, v_obtained), 0):
                        if not np.isclose(func_overlap(new_basis, b), 0):
                            # print(" already obtained")
                            flag_new = False
                            break
                    if flag_new:
                        # print("***new vector***")
                        bases_new.append(new_basis)

            # if len(bases_new) == len(bases) - 1:
            #     # TODO: the last basis is automatically determined
            #     pass
            if len(bases_new) == len(bases):
                # all bases have been determined
                break

        try:
            assert len(bases_new) == len(bases)
        except AssertionError:
            # TODO: warning
            warnings.warn("The bases not updated")
            return None
        else:
            return bases_new

    def _recombine_cf_states(self):
        self._print_verbose("\n_recombine_cf_states")

        def func_matrix_element(_vec_i, _vec_j, _op):
            return _vec_i.conj().T @ _op @ _vec_j

        def func_overlap(_vec1, _vec2):
            return _vec1.conj().T @ _vec2

        for cf_index, cf in enumerate(self._cf_states):
            # print(cf.name, cf.dim)
            if cf.dim == 1:
                continue

            # matrix elements of the symmetry operator
            vec = [cf.vecs[:, i] for i in range(cf.dim)]  # list of column vector
            vec_new = self._recombine_bases_core(vec, func_matrix_element, func_overlap)
            assert isinstance(vec_new, list)
            assert len(vec_new) == cf.dim

            # update
            if vec_new is not None:
                # vecs_new = np.hstack(vec_new)
                vecs_new = np.array(vec_new).reshape((cf.dim, -1)).T  # convert to column vectors
                assert vecs_new.shape == (self.dim, cf.dim)
                # print(vecs_new)
                self._cf_states[cf_index] = PointGroup.EigenStates(cf.name, cf.dim, vecs_new)

    def _recombine_multipoles(self):
        self._print_verbose("\n_recombine_multipoles")

        def func_matrix_element(_mat_i, _mat_j, _op):
            return matrix_overlap(_mat_i, _op @ _mat_j @ _op.conj().T)

        def func_overlap(_mat1, _mat2):
            return matrix_overlap(_mat1, _mat2)

        for irp, m_irp in list(self._multipoles.items()):
            self._print_verbose(irp)
            # for m in m_irp:
            #     print(m)
            c = Counter([m.rank for m in m_irp])
            self._print_verbose(c)

            for rank, count in list(c.items()):
                if count == 1:
                    continue
                self._print_verbose("rank =", rank)
                mat_list = [m.mat for m in m_irp if m.rank == rank]
                assert len(mat_list) == count

                mat_list_new = self._recombine_bases_core(mat_list, func_matrix_element, func_overlap)

                # update
                if mat_list_new is not None:
                    # get index of operators with a given rank
                    indices = [i for i, m in enumerate(m_irp) if m.rank == rank]
                    for i, mat in zip(indices, mat_list_new):
                        m_irp[i].mat = mat
