from linicrypt_solver.field import GF
import numpy as np
from galois import FieldArray


def stack_matrices(A: FieldArray, B: FieldArray, axis=0) -> FieldArray:
    return GF(np.concatenate((A, B), axis=axis))


def embed_left(A: FieldArray, dim: int) -> FieldArray:
    assert A.shape[1] <= dim
    zeros = GF.Zeros((A.shape[0], dim - A.shape[1]))
    return GF(np.concatenate((A, zeros), axis=1))


def embed_right(A: FieldArray, dim: int) -> FieldArray:
    assert A.shape[1] <= dim
    zeros = GF.Zeros((A.shape[0], dim - A.shape[1]))
    return GF(np.concatenate((zeros, A), axis=1))
