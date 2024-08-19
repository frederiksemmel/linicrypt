from dataclasses import dataclass
from linicrypt_solver.field import GF


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
        self.C = self.compute_C()

    def __str__(self):
        a, b, c, d, e, f = self.a, self.b, self.c, self.d, self.e, self.f
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

    def constraint(self):
        x = self.compute_x()
        k = self.compute_k()
        y = self.compute_y()
        return (x, k, y)

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
        match (self.a, self.b):
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
        match (self.c, self.d):
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
        match (self.e, self.f):
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


def row_with_one(d: int, n: int) -> Matrix:
    x = Matrix.zeros(1, d)
    x[0, n] = 1
    return x


# Warning: This is in the nicer merkle damgard basis per default
# TODO properly handle the basis changes as in the compression function
class MerkleDamgard:
    def __init__(self, f: PGVComporessionFunction, n: int):
        f.to_output_1()
        self.f = f
        self.n = n
        self.dimension = n * 2 + 1
        self.B = Identity(n * 2 + 1)

    @property
    def I(self):  # noqa: E743
        h0 = row_with_one(self.dimension, 0)
        rows = [h0]
        for i in range(1, self.dimension, 2):
            mi = row_with_one(self.dimension, i)
            rows.append(mi)
        return Matrix(rows)

    @property
    def O(self):  # noqa: E743
        o1 = row_with_one(self.dimension, 0)
        on = row_with_one(self.dimension, self.dimension - 1)
        return Matrix([o1, on])

    def c(self, i: int):  # noqa: E743
        assert i <= self.n
        assert i > 0
        k = Matrix.zeros(1, self.dimension)
        k[0, (i - 1) * 2 : i * 2 + 1] = self.f.k
        x = Matrix.zeros(1, self.dimension)
        x[0, (i - 1) * 2 : i * 2 + 1] = self.f.x
        y = Matrix.zeros(1, self.dimension)
        y[0, (i - 1) * 2 : i * 2 + 1] = self.f.y
        return (x, k, y)

    def permute_eqs(self, i: int, j: int, B: Matrix):
        xi, ki, yi = self.c(i)
        xj, kj, yj = self.c(j)
        eqs = [
            Eq(xi * B, xj),
            Eq(xj * B, xi),
            Eq(ki * B, kj),
            Eq(kj * B, ki),
            Eq(yi * B, yj),
            Eq(yj * B, yi),
        ]
        return eqs

    def collapse_eqs(self, i: int, j: int, B: Matrix):
        xi, ki, yi = self.c(i)
        xj, kj, yj = self.c(j)
        eqs = [
            Eq(xi * B, xi),
            Eq(xj * B, xi),
            Eq(ki * B, ki),
            Eq(kj * B, ki),
            Eq(yi * B, yi),
            Eq(yj * B, yi),
        ]
        return eqs

    def cycle_eqs(self, B: Matrix):
        x1, k1, y1 = self.c(1)
        xn, kn, yn = self.c(self.n)
        eqs = [
            Eq(xn * B, x1),
            Eq(kn * B, k1),
            Eq(yn * B, y1),
        ]
        for i in range(1, self.n):
            xi, ki, yi = self.c(i)
            xi_1, ki_1, yi_1 = self.c(i + 1)
            eqs += [
                Eq(xi * B, xi_1),
                Eq(ki * B, ki_1),
                Eq(yi * B, yi_1),
            ]

        return eqs

    # This doesnt make sense yet, but there has to be something there...
    # def flip_permute(self, i: int, j: int, B: Matrix):
    #     xi, ki, yi = self.c(i)
    #     xj, kj, yj = self.c(j)
    #     eqs = [
    #         Eq(xi * B, xj),
    #         Eq(xj * B, xi),
    #         Eq(ki * B, kj),
    #         Eq(kj * B, ki),
    #         Eq(yi * B, yj),
    #         Eq(yj * B, yi),
    #     ]
    #     return eqs

    def output_invariant(self, B: Matrix):
        return Eq(self.O @ B, self.O)

    def output_h0(self, B: Matrix):
        o1 = row_with_one(self.dimension, 0)
        on = row_with_one(self.dimension, self.dimension - 1)
        return Eq(on @ B, o1)

    def iv_invariant(self, B: Matrix):
        o1 = row_with_one(self.dimension, 0)
        return Eq(o1 @ B, o1)
