from dataclasses import dataclass

from loguru import logger

from linicrypt_solver.field import GF
from linicrypt_solver.ideal_cipher import ConstraintE
from linicrypt_solver.solvable import Constraints


@dataclass
class PGVParams:
    a: int
    b: int
    c: int
    d: int
    e: int
    f: int


class PGVComporessionFunction:
    def __init__(self, params: PGVParams):
        self.params = params
        self.I = self.compute_I()
        self.O = self.compute_O()
        self.C = self.construct_constraint()

    def __str__(self):
        a, b, c, d, e, f = (
            self.params.a,
            self.params.b,
            self.params.c,
            self.params.d,
            self.params.e,
            self.params.f,
        )
        return f"E({c}h + {d}m, {e}h + {f}m) + {a}h + {b}m\t PGV: {self.pgv_category()}, BRS: {self.brs_category()}"

    def compute_I(
        self,
    ):  # noqa: E743
        return GF([[1, 0, 0], [0, 1, 0]])

    def compute_O(self):  # noqa: E743
        return GF([[self.params.a, self.params.b, 1]])

    def compute_x(self):  # noqa: E743
        return GF([[self.params.c, self.params.d, 0]])

    def compute_k(self):  # noqa: E743
        return GF([[self.params.e, self.params.f, 0]])

    def compute_y(self):  # noqa: E743
        return GF([[0, 0, 1]])

    def construct_constraint(self):
        x = self.compute_x()
        k = self.compute_k()
        y = self.compute_y()
        return ConstraintE(x, k, y)

    def linicrypt_is_secure(self):
        a, b, c, d, e, f = (
            self.params.a,
            self.params.b,
            self.params.c,
            self.params.d,
            self.params.e,
            self.params.f,
        )
        return (
            (0, 0) not in {(a, b), (c, d), (e, f)}
            and (a, b) != (c, d)
            and (c, d) != (e, f)
        )

    def pgv_choice_of_ff(self):
        match (self.params.a, self.params.b):
            case (0, 0):
                return "V", 0
            case (0, 1):
                return "X_i", 1
            case (1, 0):
                return "H_{i-1}", 2
            case (1, 1):
                return "X_i + H_{i-1}", 3
        raise ValueError

    def pgv_choice_of_k(self):
        match (self.params.c, self.params.d):
            case (0, 0):
                return "V", 3
            case (0, 1):
                return "X_i", 0
            case (1, 0):
                return "H_{i-1}", 1
            case (1, 1):
                return "X_i + H_{i-1}", 2
        raise ValueError

    def pgv_choice_of_p(self):
        match (self.params.e, self.params.f):
            case (0, 0):
                return "V", 3
            case (0, 1):
                return "X_i", 0
            case (1, 0):
                return "H_{i-1}", 1
            case (1, 1):
                return "X_i + H_{i-1}", 2
        raise ValueError

    def pgv_category(self):
        _, i_ff = self.pgv_choice_of_ff()
        _, i_k = self.pgv_choice_of_k()
        _, i_p = self.pgv_choice_of_p()

        submatrices = [
            [
                [("-", None), ("B", 13), ("B", 25), ("-", None)],
                [("D", 1), ("-", None), ("D", 26), ("-", None)],
                [("B", 2), ("B", 14), ("F", 27), ("F", 1)],
                [("-", None), ("-", None), ("D", 28), ("-", None)],
            ],
            [
                [("-", None), ("B", 15), ("B", 29), ("-", None)],
                [("✓", 3), ("D", 16), ("✓", 30), ("D", 12)],
                [("FP", 4), ("FP", 17), ("B", 31), ("B", 43)],
                [("-", None), ("D", 18), ("B", 32), ("-", None)],
            ],
            [
                [("P", 5), ("FP", 19), ("FP", 33), ("P", 44)],
                [("D", 6), ("-", None), ("D", 34), ("-", None)],
                [("FP", 7), ("FP", 20), ("B", 35), ("B", 45)],
                [("D", 8), ("-", None), ("D", 36), ("-", None)],
            ],
            [
                [("P", 9), ("FP", 21), ("FP", 37), ("P", 46)],
                [("✓", 10), ("D", 22), ("✓", 38), ("D", 47)],
                [("B", 11), ("B", 23), ("F", 39), ("F", 48)],
                [("P", 12), ("D", 24), ("F", 40), ("D", 49)],
            ],
        ]

        # print(i_ff, i_k, i_p)
        pgv_category, _ = submatrices[i_ff][i_k][i_p]
        pgv_index = 16 * i_p + 4 * i_ff + i_k + 1
        assert pgv_index <= 64
        assert pgv_index > 0
        return pgv_category, pgv_index

    def brs_category(self):
        brs_categories = (
            "abcaadeafbebfdcfcacacbebeaeaebcbcbgbcdggebgbedggaagaabgafagafbgb"
        )
        _, pgv_index = self.pgv_category()
        return brs_categories[pgv_index - 1]

    def construct_MD(self, n: int):
        constraints = Constraints([])
        for i in range(1, n + 1):
            dim = constraints.dim()
            c = self.C.embed_right(dim + 3)
            constraints = constraints.embed_left(dim + 3)
            # TODO collapse the first input of this compression function with the output of the ent constraints
            constraints.add(c)
        logger.debug(f"constraints are {constraints}")
        return constraints
