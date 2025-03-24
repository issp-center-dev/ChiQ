import numpy as np
import math
from itertools import product
from fractions import Fraction
from sympy.physics.quantum.cg import CG
from sympy import S

def commutator(mat1, mat2):
    return mat1*mat2 - mat2*mat1

def anticommutator(mat1, mat2):
    return mat1*mat2 + mat2*mat1

def is_commute(mat1, mat2):
    return np.allclose(mat1*mat2, mat2*mat1)

def is_anticommute(mat1, mat2):
    return np.allclose(mat1*mat2, -mat2*mat1)

def is_hermitian(mat):
    return np.allclose(mat, np.conjugate(mat.T))
    # dif = np.linalg.norm( mat - np.conjugate(mat.T) )  # norm of A-A^h
    # if dif < 1e-10:
    #     return True
    # else:
    #     return False

def is_unitary(mat):
    # assert type(mat)
    return np.allclose( np.matrix(mat)*np.matrix(mat.T), np.identity(mat.shape[0]) )

def eigen(A):
    eigenValues, eigenVectors = np.linalg.eig(A)
    idx = eigenValues.argsort()[::-1]
    eigenValues = eigenValues[idx]
    eigenVectors = eigenVectors[:,idx]
    return eigenValues, eigenVectors

def matrix_prop_fac(_A, _B):
    """Find a such that a*A=B.

    Return a, residual error
    """

    # transform matrix to ndarray
    A = np.array(_A)
    B = np.array(_B)

    AA = A * A  # element-wise product
    AB = A * B
    a = np.sum(AB) / np.sum(AA)
    temp = a*A-B
    res = np.sum(temp*temp)
    return a, res

def find_coef(_As, _B):
    """
    Find {a_i} that fulfills sum_i a_i*A_i = B

    Input
        As: (list of np.array) Matrices
        B : (np.array) Matrix
    Return
        a_i: (np.array) coefficients
        res: (float)  residual error
    """
    assert isinstance(_As,list)

    # transform matrix to ndarray, since the behaviors of product are different between matrix and ndarray
    As = [ np.array(A) for A in _As ]
    B = np.array(_B)

    n=len(As)
    alpha = np.zeros((n,n))
    beta = np.zeros(n)
    gamma = np.sum(B*B) # take element-wise product, and sum up them
    for i,j in product(list(range(n)),list(range(n))):
        alpha[i,j] = np.sum(As[i]*As[j])
    for i in range(n):
        beta[i] = np.sum(As[i]*B)
    # print alpha
    # print beta
    # print gamma

    coef = np.linalg.solve(alpha,beta)
    res = np.dot(coef,np.dot(alpha,coef)) - 2.*np.dot(coef,beta) + gamma
    # print coef
    # print res
    return coef, res


class CrystalField:
    """Atomic states in crystal field potential."""

    def __init__(self, J, point_group):
        """
        J: (Integer) angular momentum
        point_group: (String) point group name
        """

        assert isinstance(J, int)
        assert point_group in ['Oh',]

        self.j = J
        self.point_group = point_group
        self.V_CF = None
        # self.V_SO = None
        self.__Vmat = None
        self.flag_SO = False
        self.flag_diag = False  # True if diagonalization has been done

        # operator Jz
        self.Jz = np.matrix(np.diag([ jz for jz in range(-J,J+1) ]))
        # print type(self.Jz)
        # print self.Jz

        # operator J+
        self.Jp = np.matrix(np.zeros( (2*J+1, 2*J+1) ))
        for jz in range(-J,J):
            i = jz+J
            self.Jp[i+1,i] = math.sqrt( (J-jz)*(J+jz+1) )
        # print type(self.Jp)
        # print self.Jp

        # operator J-
        self.Jm = np.matrix(self.Jp).T
        # print type(self.Jm)
        # print self.Jm

        # identity matrix
        self.I = np.matrix(np.identity(2*J+1))

        def commutator(mat1, mat2):
            return mat1*mat2 - mat2*mat1

        # check [Jz,J+] = J+
        assert np.allclose( commutator(self.Jz, self.Jp), self.Jp )
        # check [Jz,J-] = J-
        assert np.allclose( commutator(self.Jz, self.Jm), -self.Jm )
        # check [J+,J-] = 2 Jz
        assert np.allclose( commutator(self.Jp, self.Jm), 2*self.Jz )

        self.jsq = J*(J+1)

    def set_CF(self, B_params=None, Wx_params=None):
        self.flag_diag = False  # reset flag
        self.__Vmat = None

        if self.point_group == 'Oh':
            if B_params is not None:
                B4 = B_params['B4']
                B6 = B_params['B6']
            elif Wx_params is not None:
                W = Wx_params['W']
                x = Wx_params['x']
                #      J=0     1     2     3     4     5     6     7
                F4 = (None, None,   12,   15,   60,   14,   60,   60)
                F6 = (None, None, None,  180, 1260, 1260, 7560, 3780)
                B4 = W*x / F4[self.j]     if F4[self.j] is not None else 0
                B6 = W*(1-x) / F6[self.j] if F6[self.j] is not None else 0
            else:
                raise Exception("Either B_params or Wx_params must be given")

            op4 = self.__op40() + 5*self.__op44()
            op6 = self.__op60() - 21*self.__op64()
            self.V_CF = B4 * op4 + B6 * op6
        else:
            raise Exception("point_group %s not implemented" %point_group)

    def __op40(self):
        return 35*self.Jz**4 - (30*self.jsq-25) * self.Jz**2 + ( -6*self.jsq + 3*self.jsq**2 ) * self.I

    def __op44(self):
        return (self.Jp**4 + self.Jm**4)/2

    def __op60(self):
        return 231*self.Jz**6 - 105*(3*self.jsq-7) * self.Jz**4 + (105*self.jsq**2-525*self.jsq+294)*self.Jz**2 + ( -5*self.jsq**3 + 40*self.jsq**2 - 60*self.jsq ) * self.I

    def __op64(self):
        return anticommutator( (11*self.Jz**2 - self.jsq * self.I - 38*self.I), (self.Jp**4 + self.Jm**4) ) / 4

    def set_SO(self, _lambda, order_of_spin_components):
        """
        order_of_spin_components: 'up_dn' or 'dn_up'
        """
        self.flag_diag = False  # reset flag
        self.flag_SO = True

        dim = 2*self.j+1
        self.V_SO = np.zeros((2*dim,2*dim))
        if order_of_spin_components == 'up_dn':
            self.V_SO[0:dim, 0:dim] += self.Jz[:,:]/2.  # L_z.S_z
            self.V_SO[dim:, dim:] += -self.Jz[:,:]/2.   # L_z.S_z
            self.V_SO[0:dim, dim:] += self.Jm[:,:]/2.  # L_-.S_+ / 2
            self.V_SO[dim:, 0:dim] += self.Jp[:,:]/2.  # L_+.S_- / 2
        elif order_of_spin_components == 'dn_up':
            self.V_SO[0:dim, 0:dim] += -self.Jz[:,:]/2.  # L_z.S_z
            self.V_SO[dim:, dim:] += self.Jz[:,:]/2.   # L_z.S_z
            self.V_SO[0:dim, dim:] += self.Jp[:,:]/2.  # L_+.S_- / 2
            self.V_SO[dim:, 0:dim] += self.Jm[:,:]/2.  # L_-.S_+ / 2
        else:
            raise ValueError(order_of_spin_components)
        self.V_SO *= _lambda
        # print self.V_SO

    def unset_SO(self):
        self.flag_diag = False  # reset flag
        self.flag_SO = False
        self.V_SO = None

    def set_Vmat_external(self, _matrix):
        self.flag_diag = False  # reset flag

        assert is_hermitian(_matrix)
        self.__Vmat = np.array(_matrix)

        # check matrix size and set SO flag
        if _matrix.shape == (2*self.j+1, 2*self.j+1):
            self.flag_SO = False
            self.V_CF = self.__Vmat
        elif _matrix.shape == (2*(2*self.j+1), 2*(2*self.j+1)):
            self.flag_SO = True
            self.V_CF = None
        else:
            Exception("Matrix size incorrect:", _matrix.shape)
            return

    def __diagonalize(self):
        if self.flag_diag:  # already diagonalized
            return

        self.flag_diag = True
        Vmat = self.get_Vmat()
        self.levels, self.vecs = eigen(Vmat)

    def get_Vmat(self):
        if self.__Vmat is not None:
            return self.__Vmat.copy()
        else:
            if self.flag_SO:  # w/ SO
                Vmat = self.V_SO.copy()
                dim = 2*self.j+1
                Vmat[0:dim, 0:dim] += self.V_CF[:,:]
                Vmat[dim:, dim:] += self.V_CF[:,:]
            else:  # w/o SO (only CF potential)
                Vmat = self.V_CF.copy()
            assert is_hermitian(Vmat)
            return Vmat

    def get_levels(self):
        self.__diagonalize()
        return self.levels

    def get_vectors(self):
        self.__diagonalize()
        return self.vecs

    def print_levels(self):
        self.__diagonalize()

        print("Energy levels:")
        # print self.levels
        E0 = self.levels[-1]
        for i,x in enumerate(self.levels):
            print(" %2d : %12.6f (%12.6f)" %(self.levels.shape[0]-i-1, x, x-E0))

    def print_vectors(self, LS_or_J='LS', rationalize=False):
        """
        LS_or_J = 'LS' or 'J'
        """
        self.__diagonalize()

        print("Eigenvectors:")

        # basis
        if self.flag_SO:
            dim = 2*self.j+1
            lz = list(range(-self.j,self.j+1)) * 2
            sz = [S(-1)/2]*dim + [S(1)/2]*dim
            j = [self.j+S(1)/2]*(dim+1) + [self.j-S(1)/2]*(dim-1)
            jz = list(range(-2*self.j-1, 2*self.j+3, 2)) + list(range(-2*self.j+1, 2*self.j+1, 2))  # 2*jz
            jz = [ S(x)/2 for x in jz ]
            assert len(lz)==2*dim
            assert len(sz)==2*dim
            assert len(j)==2*dim
            assert len(jz)==2*dim
        else:
            lz = list(range(-self.j,self.j+1))

        # construct U matrix
        if LS_or_J == 'LS':
            print(" Lz =", lz)
            if self.flag_SO:
                print(" Sz =", sz)
            U = np.matrix(np.identity(self.vecs.shape[0]))

        elif LS_or_J == 'J':
            if not self.flag_SO:
                raise Exception("")
            print(" J  =", j)
            print(" Jz =", jz)
            M=len(lz)
            N=len(j)
            U = np.matrix(np.zeros((M,N))) # < J | LS >
            for m in range(M):
                for n in range(N):
                    # < L Lz S Sz | J Jz >
                    cg = CG(self.j, lz[m], S(1)/2, sz[m], j[n], jz[n])
                    U[m,n] = cg.doit().evalf()
            U = U.H

        def float_to_frac(x):
            # print x,
            tol = 1e-10
            if abs(x) < tol:
                return "0,"
            elif abs(x-1) < tol:
                return "1,"
            elif abs(x+1) < tol:
                return "-1,"
            else:
                frac = Fraction(x**2).limit_denominator()
                r = "sqrt(%s)," %frac
                if x>0:
                    return r
                else:
                    return "-%s" %r

        def float_to_str(x):
            tol = 1e-10
            if abs(x) < tol:
                return " 0   "
            elif abs(x-1) < tol:
                return " 1   "
            elif abs(x+1) < tol:
                return "-1   "
            else:
                return "%5.2lf" %x

        if rationalize==True:
            str_x = float_to_frac
        else:
            str_x = float_to_str

        # array = np.array(self.vecs)
        array = np.array(U*self.vecs)
        for i in range(array.shape[1]):
            vec = array[:,i]
            # print vec
            print(" %2d :  (" %(self.levels.shape[0]-i-1), end=' ')
            for x in vec:
                # print x,
                print("%s" %str_x(x), end=' ')
            print(")")

    def projection(self):
        """Classify eigenvectors into irr reps by projection"""
        self.__diagonalize()

        print("Projection to CEF basis")
        def close_basis(_array):
            # array = np.array(_array)[0]  # matrix to ndarray
            # print array.shape
            index = np.argmax(_array)
            # print array[index]
            if _array[index] > 0.999:
                return index
            else:
                return None

        def indexsort(_array):
            idx = _array.argsort()[::-1]
            return idx

        proj = self.vecs.T * self.__U()
        proj = np.array(proj)  # matrix to ndarray
        proj = np.array([ [x**2 for x in row] for row in proj ])  # x**2 for each element x
        for i in range(proj.shape[1]):
            print(" %2d :" %(self.levels.shape[0]-i-1), end=' ')
            weights = proj[i,:]

            # j = close_basis(weights)
            # if j is not None:
            #     print self.basis_names[j]
            # else:
            #     print weights

            # sort by weights
            idx = indexsort(weights)
            weights_sorted = weights[idx]
            names = np.array(self.basis_names)[idx]
            # print weights_sorted
            for i,w in enumerate(weights_sorted):
                if w>1e-3:
                    print(" %s (%.3lf)," %(names[i], w), end=' ')
            print("")

    def __U(self):
        """Unitary matrix which diagonalizes V_CF (set of column vectors)"""

        if self.point_group == 'Oh':
            if self.j == 0:
                self.basis_names = ['s']
                ut = [[1.],
                     ]
            elif self.j == 1:
                self.basis_names = ['px', 'py', 'pz']
                ut = [[1./2., 0, 1./2.,],  # px
                      [1./2., 0, -1./2.],  # py
                      [0, 1., 0],  # pz
                     ]
            elif self.j == 2:
                self.basis_names = ['eg', 'eg', 't2g', 't2g', 't2g']
                # self.basis_names = ['Gamma3', 'Gamma3', 'Gamma5', 'Gamma5', 'Gamma5']
                ut = [[0, 0, 1., 0, 0],  # Gamma3
                      [1./2., 0, 0, 0, 1./2.],  # Gamma3
                      [0, 0, 0, 1., 0],  # Gamma5
                      [0, 1., 0, 0, 0],  # Gamma5
                      [-1./2., 0, 0, 0, 1./2.],  # Gamma5
                     ]
            elif self.j == 3:
                self.basis_names = ['Gamma2', 'Gamma4', 'Gamma4', 'Gamma4', 'Gamma5', 'Gamma5', 'Gamma5']
                ut = [[0, -1./2., 0, 0, 0, 1./2., 0],  # Gamma2
                      [0, 0, 3./8., 0, 0, 0, 5./8.],  # Gamma4
                      [5./8., 0, 0, 0, 3./8., 0, 0],  # Gamma4
                      [0, 0, 0, 1., 0, 0, 0],  # Gamma4
                      [0, 1./2., 0, 0, 0, 1./2., 0],  # Gamma5
                      [0, 0, 5./8., 0, 0, 0, -3./8.],  # Gamma5
                      [3./8., 0, 0, 0, -5./8., 0, 0],  # Gamma5
                     ]
            elif self.j == 4:
                self.basis_names = ['Gamma1', 'Gamma3', 'Gamma3', 'Gamma4', 'Gamma4', 'Gamma4', 'Gamma5', 'Gamma5', 'Gamma5']
                ut = [[5./24., 0, 0, 0, 7./12., 0, 0, 0, 5./24.],  # Gamma1
                      [0, 0, 1./2., 0, 0, 0, 1./2., 0, 0],  # Gamma3
                      [7./24., 0, 0, 0, -5./12., 0, 0, 0, 7./24.],  # Gamma3
                      [0, 1./8., 0, 0, 0, 7./8., 0, 0, 0],  # Gamma4
                      [0, 0, 0, 7./8., 0, 0, 0, 1./8., 0],  # Gamma4
                      [-1./2., 0, 0, 0, 0, 0, 0, 0, 1./2.],  # Gamma4
                      [0, 0, 0, 1./8., 0, 0, 0, -7./8., 0],  # Gamma5
                      [0, -7./8., 0, 0, 0, 1./8., 0, 0, 0],  # Gamma5
                      [0, 0, 1./2., 0, 0, 0, -1./2., 0, 0],  # Gamma5
                     ]
            else:
                raise Exception("J=%d not implemented" %self.j)
        else:
            raise Exception("point_group %s not implemented" %point_group)

        assert len(self.basis_names) == 2*self.j+1
        ut_sgn = np.sign(ut)
        ut = [ [math.sqrt(abs(x)) for x in row] for row in ut ]  # sqrt(|x|) for each element x
        ut = ut*ut_sgn  # element-wise product
        assert is_unitary(np.array(ut))
        return np.matrix(ut).T


    def extract_SO(self):
        assert self.flag_SO == True

        print("Extract SOC constant, and decompose Vmat into V_CF and V_SO")

        Vmat = self.get_Vmat()
        dim = 2*self.j+1
        M11 = Vmat[0:dim, 0:dim]
        M22 = Vmat[dim:, dim:]
        M12 = Vmat[0:dim, dim:]
        # M21 = Vmat[dim:, 0:dim]

        # so_z, res_z = matrix_prop_fac(self.Jz, M22-M11)
        # so_x, res_x = matrix_prop_fac(self.Jp, 2.*M12)
        [so_z,], res_z = find_coef([self.Jz,], M22-M11)
        [so_x,], res_x = find_coef([self.Jp,], 2.*M12)

        print(" SO coupling (residual error)")
        print("  z:", so_z, "(%.2e)" %res_z)
        print("  x:", so_x, "(%.2e)" %res_x)

        self.V_CF = np.array(M11 + M22)/2.
        return (so_z+so_x)/2., self.V_CF.copy()

    def extract_CF_params(self):
        if self.V_CF is None:
            Exception("V_CF is not given. Call set_Vmat_external and extract_SO (if flag_SO) before this method")

        print("Extract CF parameters")

        if self.point_group == 'Oh':
            op4 = self.__op40() + 5*self.__op44()
            op6 = self.__op60() - 21*self.__op64()

            #      J=0     1     2     3     4     5     6     7
            F4 = (None, None,   12,   15,   60,   14,   60,   60)
            F6 = (None, None, None,  180, 1260, 1260, 7560, 3780)

            # find B values that fulfills
            # self.V_CF = B4 * op4 + B6 * op6
            if self.j>=3:
                [B4, B6], res = find_coef([op4, op6], self.V_CF)
                W = B4*F4[self.j] + B6*F6[self.j]
                x = B4*F4[self.j]/W
            elif self.j==2:
                B6=0
                [B4,], res = find_coef([op4,], self.V_CF)
                W = B4*F4[self.j]
                x = 1
            else:
                B4=0
                B6=0
                W=0
                x=0
            B_params = {'B4': B4, 'B6': 0}
            Wx_params = {'W': W, 'x': x}
        else:
            raise Exception("point_group %s not implemented" %point_group)

        print("", B_params)
        print("", Wx_params)
        print(" (residual error: %.2e)" %res)
        return B_params, res


if __name__ ==  "__main__":
    L=3
    CF = CrystalField(J=L, point_group='Oh')
    W=1.
    x=-0.4

    print("\n==== w/o SO")
    CF.set_CF(Wx_params={'W':W, 'x':x})
    CF.print_levels()
    CF.print_vectors(rationalize=True)
    CF.projection()

    print("\n==== w/ SO")
    CF.set_SO(100., order_of_spin_components='dn_up')
    CF.print_levels()
    # CF.print_vectors()
    CF.print_vectors(rationalize=False, LS_or_J='J')
    vmat = CF.get_Vmat()

    print("\n==== analyze Vmat")
    CF_ext = CrystalField(J=L, point_group='Oh')
    CF_ext.set_Vmat_external(vmat)
    CF_ext.extract_SO()
    CF_ext.extract_CF_params()

    # test
    # val,vec = eigen([[0,1],[1,0]])
    # print val
    # # print vec
    # print_vectors(vec)
