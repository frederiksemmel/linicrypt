import sys

from loguru import logger
from itertools import product

from linicrypt_solver.algebraic_representation import AlgebraicRep
from linicrypt_solver.field import GF
from linicrypt_solver.merkle_damgard import PGVComporessionFunction, PGVParams
from linicrypt_solver.solvable import Constraints

# Configure logger to print the log message on a new line
logger.remove()  # Remove the default handler
logger.add(
    sink=sys.stderr,
    level="INFO",
    format="{time:YYYY-MM-DD HH:mm:ss} | {level}\n{message}",  # Add newline before {message}
)


def test_cr():
    constraints = Constraints.from_repr(
        [
            ([1, 0, 0, 0, 0], [0, 0, 1, 0, 0]),
            ([0, 0, 1, 0, 0], [0, 0, 0, 1, 0]),
            ([0, 1, 0, 0, 0], [0, 0, 0, 0, 1]),
        ]
    )
    fixing = GF([[1, 0, 0, 0, 0], [0, 1, 0, 0, 0]])
    output = GF([[0, 0, 0, 1, 1]])
    program = AlgebraicRep(constraints, fixing, output)
    # print(program)
    # program.cs.is_solvable(fixing=program.fixing)
    print(f"The program\n{program}\n")
    print(f"is CR: {program.is_collision_resistant()}")
    print(f"is 2PR: {program.is_second_preimage_resistant()}")


def test_MD_with(a, b, c, d, e, f):
    params = PGVParams(a, b, c, d, e, f)
    pgv_f = PGVComporessionFunction(params)
    # if pgv_f.pgv_category()[0] != "B":
    #     continue
    print(pgv_f)
    n = 3
    H_n = pgv_f.construct_MD(n)
    logger.debug(f"H_n:\n{H_n}")
    # print(H_n.cs.is_solvable(fixing=H_n.fixing))
    print(f"Collision resistant: {H_n.is_collision_resistant()}")
    # print(f"Collision resistant: {H_n.is_second_preimage_resistant()}")
    # break


def test_MD():
    for params in product([0, 1], repeat=6):
        test_MD_with(*params)


if __name__ == "__main__":
    test_cr()
    # test_MD()
    # test_MD_with(1, 0, 1, 1, 1, 0)
