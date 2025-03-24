from chiq.matrix_dict import MatrixDict
import numpy as np


def make_hermitian(matrix):
    assert matrix.shape[0] == matrix.shape[1]
    return (matrix + matrix.conj().T)/2.


def test_matrix_dict():
    # set up matrix
    sizes = [3,4,5,6]
    A = {}

    def set_matrix(i,j):
        A[(i,j)] = np.random.randn(sizes[i], sizes[j])

    set_matrix(0,0)
    set_matrix(1,1)
    set_matrix(0,2)
    # set_matrix(2,0)
    set_matrix(2,2)
    set_matrix(3,3)
    set_matrix(0,3)
    set_matrix(3,0)

    A[(1,1)] = make_hermitian(A[(1,1)])

    # diagonalize A
    BM = MatrixDict(A, verbose=True)
    eigen = BM.get_eigen()

    # compute eigen by U^{-1}*A*U and check it
    assert np.allclose(eigen, BM.compute_eigen(A)[0])

    # reconstruct the matrix from eigen
    B = BM.transform_diag(eigen)

    # compare A and B
    print("A.keys =", list(A.keys()))
    print("B.keys =", list(B.keys()))
    assert set(A.keys()) == set(B.keys())
    assert all( A[key].shape == B[key].shape for key in list(A.keys()) )
    assert all( np.allclose(A[key], B[key]) for key in list(A.keys()) )

    # trivial check
    assert BM.check_if_diagonalized(A)  # A is diagonalized by U
    set_matrix(0,0)  # change A
    assert BM.check_if_diagonalized(A) == False  # A is NOT diagonalized by U
