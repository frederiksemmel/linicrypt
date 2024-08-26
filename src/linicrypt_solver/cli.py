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
    level="WARNING",
    format="{level}\n{message}",  # Add newline before {message}
)


def example_no_nonces() -> AlgebraicRep:
    constraints = Constraints.from_repr(
        [
            ([1, 0, 0, 0, 0], [0, 0, 1, 0, 0]),
            ([0, 0, 1, 0, 0], [0, 0, 0, 1, 0]),
            ([0, 1, 0, 0, 0], [0, 0, 0, 0, 1]),
        ]
    )
    fixing = GF([[1, 0, 0, 0, 0], [0, 1, 0, 0, 0]])
    output = GF([[0, 0, 0, 1, 1]])
    return AlgebraicRep(constraints, fixing, output)


def running_example() -> AlgebraicRep:
    constraints = Constraints.from_repr(
        [
            ([1, 0, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0]),
            ([0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 1, 0]),
            ([0, 1, 1, 1, 0, 0], [0, 0, 0, 0, 0, 1]),
        ]
    )
    fixing = GF([[1, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0]])
    output = GF([[1, 0, 0, 0, 0, 1], [0, 0, 0, 0, 1, 0]])
    return AlgebraicRep(constraints, fixing, output)


def test_cr(program: AlgebraicRep):
    attacks = list(program.all_maximal_collision_attacks())
    if len(attacks) == 0:
        print(f"The program\n{program}\nis Collision Resistant")
    else:
        print(f"The program\n{program}\nhas the following maximal attacks")
        for attack in attacks:
            print(f"{attack}")
    # print(f"2PR attacks:\n{list(program.list_second_preimage_attacks())}")


def test_MD_with(a, b, c, d, e, f):
    params = PGVParams(a, b, c, d, e, f)
    pgv_f = PGVComporessionFunction(params)
    # if pgv_f.pgv_category()[0] != "B":
    #     continue
    print(pgv_f)
    n = 2
    H_n = pgv_f.construct_MD(n)
    logger.debug(f"H_n:\n{H_n}")
    test_cr(H_n)


def test_MD():
    for params in product([0, 1], repeat=6):
        test_MD_with(*params)


if __name__ == "__main__":
    # test_cr(running_example())
    test_MD()
    # test_MD_with(1, 0, 1, 1, 1, 0)
