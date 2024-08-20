import sys

from loguru import logger

from linicrypt_solver.field import GF
from linicrypt_solver.merkle_damgard import PGVComporessionFunction, PGVParams
from linicrypt_solver.solvable import Constraints
from linicrypt_solver.utils import embed_left, embed_right

# Configure logger to print the log message on a new line
logger.remove()  # Remove the default handler
logger.add(
    sink=sys.stderr,
    level="DEBUG",
    format="{time:YYYY-MM-DD HH:mm:ss} | {level}\n{message}",  # Add newline before {message}
)


def test_cr():
    C = Constraints.from_repr(
        [
            ([1, 0, 0, 0, 0], [0, 0, 1, 0, 0]),
            ([0, 0, 1, 0, 0], [0, 0, 0, 1, 0]),
            ([0, 1, 0, 0, 0], [0, 0, 0, 0, 1]),
        ]
    )
    output = GF([[0, 0, 0, 2, 2]])
    S = GF(
        [
            [1, 0, 0, 0, 0],
            [0, 1, 0, 0, 0],
            [0, 0, 1, 0, 0],
            [0, 0, 0, 1, 0],
            [0, 0, 0, 0, 1],
            [1, 0, 0, 0, 0],
            [0, 1, 0, 0, 0],
            [0, 0, 1, 0, 0],
            [0, 0, 0, 1, 0],
            [0, 0, 0, 0, 1],
        ]
    )
    C_join = C.construct_joined()
    dim = C_join.dim()
    output_collapse = embed_left(output, dim) - embed_right(output, dim)
    f = output_collapse.null_space().transpose()

    # Robust way to compute the preimage of S
    # Annihlator of S called S^0 are the dual vectors that are zero on S
    S_0 = S.left_null_space()
    # This is f^*(S^0)
    f_pullback_S_0 = S_0 @ f
    # We have f^-1(S) = f^*(S^0))^0. Because S is in the image of f, this f(f^-1(S)) = S
    preimage_S = f_pullback_S_0.null_space().transpose()
    assert (f @ preimage_S == S).all()

    C_joined_f = C_join.map(f)
    subspaces = C_joined_f.find_solvable_subspaces_outside(preimage_S)
    for part, subspace in subspaces:
        print(part)
        print(subspace)


def test_MD():
    params = PGVParams(1, 1, 0, 1, 1, 0)
    pgv_f = PGVComporessionFunction(params)
    n = 3
    pgv_f.construct_MD(n)


if __name__ == "__main__":
    # test_cr()
    test_MD()
