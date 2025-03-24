

import numpy as np
from itertools import product
import math
import sys

zero_tol=1e-8

def is_hermitian(matrix):
    return np.allclose(matrix, matrix.conj().T)

def is_zeromatrix(matrix):
    return True if np.linalg.norm(matrix) < zero_tol else False

def is_diagonal(matrix):
    for i,j in product(list(range(matrix.shape[0])),list(range(matrix.shape[1]))):
        if i != j:
            if abs(matrix[i,j]) > zero_tol:
                return False
    return True

def norm_offdiagonal(matrix):
    matrix_offd = matrix.copy()
    for i in range(min(matrix.shape[0],matrix.shape[1])):
        matrix_offd[i,i] = 0
    return np.linalg.norm(matrix_offd)

def max_offdiagonal(matrix):
    matrix_offd = np.absolute(matrix)  # element-wise abs
    for i in range(min(matrix.shape[0],matrix.shape[1])):
        matrix_offd[i,i] = 0
    return matrix_offd.max()

def diagonalize(matrix):
    assert isinstance(matrix, np.ndarray)
    assert matrix.ndim == 2
    assert matrix.shape[0] == matrix.shape[1]  # square matrix

    if is_hermitian(matrix):
        # print "*** hermitian"
        eig, u = np.linalg.eigh(matrix)
        uinv = u.conj().T
    else:
        # print "*** non-hermitian"
        eig, u = np.linalg.eig(matrix)
        uinv = np.linalg.inv(u)

    return u, eig, uinv

def str_vec_format(v):
    """
    return string for printing vector in a readable format
    """

    r = ""
    for z in v:
        # str += "  %6.3lf" % z.real if abs(z.real) > 1e-4 else "   0    "
        # str += " %+6.3lfj" % z.imag if abs(z.imag) > 1e-4 else "        "
        r += "  %4.1lf" % z.real if abs(z.real) > 1e-4 else "   0  "
        r += " %+4.1lfj" % z.imag if abs(z.imag) > 1e-4 else "      "
    return r

def str_vec_format_2(v):
    """
    return string for printing vector as a matrix
    """

    # get matrix_size = sqrt(vector_size)
    size = v.shape[0]
    size_sqrt = math.sqrt(size)
    assert size_sqrt.is_integer()
    n = int(size_sqrt)
    # print n

    mat = v.reshape((n,n))
    r = ""
    for i in range(n):
        for j in range(n):
            z = mat[i,j]
            r += "  %6.3lf" % z.real if abs(z.real) > 1e-4 else "   0    "
            r += " %+6.3lfj" % z.imag if abs(z.imag) > 1e-4 else "        "
            # str += "  %4.1lf" % z.real if abs(z.real) > 1e-4 else "   0  "
            # str += " %+4.1lfj" % z.imag if abs(z.imag) > 1e-4 else "      "
        r += "\n"
    return r


class MatrixDict:
    def __init__(self, _A, verbose=True):
        assert isinstance(_A, dict)
        self.verbose = verbose

        self.__find_number_of_blocks(_A)
        self.__analyze_block_structure(_A)
        self.__diagonalize(_A)

    def get_eigen(self):
        eigen = np.hstack([ self.__eig[blocks] for blocks in self.__block_structure ])
        assert eigen.shape[0] == self.__N
        return eigen

    def compute_eigen(self, _B):
        """
        compute eigenvalues of matmrix _B assuming that _B is diagonalized by U

        return:
            (np.array)eigenvalues if diagonalized, or
            None if not diagonalized
        """
        # check if block_structure is compatible
        self.__check_block_structure(_B)

        eigen_list = []
        error = 0
        for blocks in self.__block_structure:
            matrix = self.__make_matrix(_B, blocks)
            u = self.__U[blocks]
            uinv = self.__Uinv[blocks]
            # diag_matrix = np.dot( np.dot(uinv, matrix), u )
            # if is_diagonal(diag_matrix):
            #     eigen_list.append( diag_matrix.diagonal() )
            # else:
            #     return None
            uinv_matrix_u = np.dot( np.dot(uinv, matrix), u )
            eigen_list.append( uinv_matrix_u.diagonal() )
            # error = error + norm_offdiagonal(uinv_matrix_u)
            error = max(error, max_offdiagonal(uinv_matrix_u))

        eigen = np.hstack(eigen_list)
        assert eigen.shape[0] == self.__N
        # return eigen
        return eigen, error

    # def check_if_diagonalized(self, _B):
    #     """
    #     return True if matrix _B is diagonalized by U
    #     """
    #     # check if block_structure is compatible
    #     self.__check_block_structure(_B)
    #
    #     for blocks in self.__block_structure:
    #         matrix = self.__make_matrix(_B, blocks)
    #         u = self.__U[blocks]
    #         uinv = self.__Uinv[blocks]
    #         diag_matrix = np.dot( np.dot(uinv, matrix), u )
    #         if not is_diagonal(diag_matrix):
    #             return False
    #     return True

    def check_if_diagonalized(self, _B):
        """
        return True if matrix _B is diagonalized by U
        """
        eigen, error = self.compute_eigen(_B)
        return error < zero_tol

    def transform_diag(self, _diag):
        assert isinstance(_diag, np.ndarray)
        assert _diag.ndim == 1
        assert _diag.shape[0] == self.__N

        diag_dict = self.__split_into_blocks(_diag)

        B = {}
        for blocks in self.__block_structure:
            B.update( self.__transform_into_block_matrix(diag_dict[blocks], blocks) )
        return B

    def transform_matrix(self, _B):
        self.__analyze_block_structure(_B)
        # check if block_structure is compatible
        pass

    def print_eigenvectors(self, filename=None, print_as_matrix=False):
        if isinstance(filename, str):
            f = open(filename, 'w')
        else:
            f = sys.stdout

        mode = 0
        print("-----------------------------------", file=f)
        print("eigenvectors", file=f)
        for blocks in self.__block_structure:
            print("--------------", file=f)
            print("block =", blocks, file=f)
            u = self.__U[blocks]
            n_blocks = len(blocks)
            n_basis = u.shape[0]
            size_block = n_basis // n_blocks
            # print "size_block =", size_block
            for n in range(n_basis):
                print(mode, " (%d in block)" %n, file=f)
                mode += 1
                vec = u[:,n] # column vector
                if print_as_matrix:
                    for i in range(n_blocks):
                        print(str_vec_format_2(vec[size_block*i:size_block*(i+1)]), file=f)
                else:
                    print(str_vec_format(vec), file=f)
        print("---------------------------------", file=f)

        if filename is not None:
            f.close()

    def __find_number_of_blocks(self, _A):
        n=0
        for key in list(_A.keys()):
            n = max(n, key[0]+1)
            n = max(n, key[1]+1)

        self.n_blocks = n
        if self.verbose:
            print("find_number_of_blocks")
            print("  n_blocks =", self.n_blocks)

    def __analyze_block_structure(self, _A):
        block_structure = []

        for key in list(_A.keys()):
            # print key
            assert len(key) == 2
            if key[0] != key[1]:
                # if key_exist(key[0])
                key_found = False
                for blocks in block_structure:
                    if key[0] in blocks or key[1] in blocks:
                        blocks.add(key[0])
                        blocks.add(key[1])
                        # blocks.add(key)
                        key_found = True
                        break
                if not key_found:
                    block_structure.append( set(key) )

        self.__block_structure = []
        for blocks in block_structure:
            self.__block_structure.append( tuple(blocks) )
        # print self.__block_structure

        def is_block_exist(block):
            for blocks in self.__block_structure:
                if block in blocks:
                    return True
            return False

        for i in range(self.n_blocks):
            if not is_block_exist(i):
                self.__block_structure.append((i,))


        self.__block_sizes = [0] * self.n_blocks
        for key, data in list(_A.items()):
            for lr in range(2):
                if self.__block_sizes[key[lr]]:  # check consistency if exists
                    assert self.__block_sizes[key[lr]] == data.shape[lr]
                else:  # add
                    self.__block_sizes[key[lr]] = data.shape[lr]

        self.__N = np.sum(self.__block_sizes)

        if self.verbose:
            print("analyze_block_structure")
            print("  block_structure =", self.__block_structure)
            print("  block_sizes =", self.__block_sizes)
            print("  N =", self.__N)

    def __check_block_structure(self, _B):
        def get_blocks(block):
            for blocks in self.__block_structure:
                if block in blocks:
                    return blocks
            raise Exception("block '%s' not found in" %block, self.__block_structure)

        for key in list(_B.keys()):
            assert len(key) == 2
            if get_blocks(key[0]) != get_blocks(key[1]):
                raise Exception("Block structure of _B is not compatible with _A")

    def __diagonalize(self, _A):
        # self.__eig = []
        # self.__U = []
        # self.__Uinv = []
        self.__eig = {}
        self.__U = {}
        self.__Uinv = {}
        for blocks in self.__block_structure:
            matrix = self.__make_matrix(_A, blocks)
            u, eig, uinv = diagonalize(matrix)
            # self.__eig.append(eig)
            # self.__U.append(u)
            # self.__Uinv.append(uinv)
            self.__eig[blocks] = eig
            self.__U[blocks] = u
            self.__Uinv[blocks] = uinv

        # if self.verbose:
        #     print self.__eig

    def __make_matrix(self, _A, blocks):
        sizes = [ self.__block_sizes[block] for block in blocks ]
        pos = [0]
        for size in sizes:
            pos.append( pos[-1] + size )
        assert len(pos) == len(sizes)+1
        # n = sum(self.__block_sizes[block] for block in blocks)
        n = pos[-1]
        assert n == sum(self.__block_sizes[block] for block in blocks)

        if self.verbose:
            print("make_matrix")
            print("  sizes =", sizes)
            print("  pos =", pos)
            print("  n =", n)

        matrix = np.zeros((n,n), dtype=np.complex128)
        for (i1,b1), (i2,b2) in product(enumerate(blocks),enumerate(blocks)):
            # print "(i1,b1) =", i1, b1
            # print "(i2,b2) =", i2, b2
            if (b1,b2) in _A:
                # print _A[(b1,b2)]
                # print "pos[i1], pos[i2] =", pos[i1],pos[i2]
                matrix[pos[i1]:pos[i1]+sizes[i1], pos[i2]:pos[i2]+sizes[i2]] = _A[(b1,b2)][:,:]  # copy
        # print matrix
        # print matrix.shape
        return matrix

    def __split_into_blocks(self, array):
        pos = [0]
        for size in self.__block_sizes:
            pos.append( pos[-1] + size )
        assert len(pos) == len(self.__block_sizes)+1
        assert pos[-1] == self.__N

        if self.verbose:
            print("split_into_blocks")
            print("  pos =", pos)

        # print array
        cp = array.copy()
        array_dict = {}
        for blocks in self.__block_structure:
            # for block in blocks:
            #     array_cut.hstack( array[pos[block]:pos[block+1]] )
            # array_cut = np.hstack([ array2[pos[block]:pos[block+1]] for block in blocks ])
            # array_dict[blocks] = array_cut
            size = sum([ self.__block_sizes[block] for block in blocks ])
            # print size
            array_dict[blocks] = cp[0:size]
            cp = cp[size:]
        # print array_dict
        return array_dict

    def __transform_into_block_matrix(self, array, blocks):
        sizes = [ self.__block_sizes[block] for block in blocks ]
        pos = [0]
        for size in sizes:
            pos.append( pos[-1] + size )
        assert len(pos) == len(sizes)+1
        # n = sum(self.__block_sizes[block] for block in blocks)
        n = pos[-1]
        assert n == sum(self.__block_sizes[block] for block in blocks)

        u = self.__U[blocks]
        uinv = self.__Uinv[blocks]
        diag = np.diag(array)

        if self.verbose:
            print("get_block,  blocks =", blocks)

        matrix = np.dot( np.dot(u, diag), uinv )
        # print blocks
        # print matrix.shape

        B = {}
        for (i1,b1), (i2,b2) in product(enumerate(blocks),enumerate(blocks)):
            matrix_b1b2 = matrix[pos[i1]:pos[i1]+sizes[i1], pos[i2]:pos[i2]+sizes[i2]]
            if not is_zeromatrix(matrix_b1b2):
                B[(b1,b2)] = matrix_b1b2
            else:
                if self.verbose:
                    print("  zeromatrix ", b1, b2)
        return B
