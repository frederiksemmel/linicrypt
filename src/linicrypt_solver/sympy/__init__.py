from itertools import product
from dataclasses import dataclass

# from sympy import BooleanTrue
from sympy import symbols, Matrix, Eq, solve, init_printing, Identity, pprint, eye

init_printing(use_unicode=True)


@dataclass
class PGVComporessionFunction:
    a: int
    b: int
    c: int
    d: int
    e: int
    f: int
    B = Identity(3)

    def __str__(self):
        a, b, c, d, e, f = self.a, self.b, self.c, self.d, self.e, self.f
        return f"E({c}h + {d}m, {e}h + {f}m) + {a}h + {b}m\t PGV: {self.pgv_category()}, BRS: {self.brs_category()}"

    def to_canonical(self):
        self.B = Identity(3)

    def to_output_1(self):
        self.B = Matrix([[1, 0, 0], [0, 1, 0], [-self.a, -self.b, 1]])

    @property
    def I(self):  # noqa: E743
        return Matrix([[1, 0, 0], [0, 1, 0]]) @ self.B

    @property
    def O(self):  # noqa: E743
        return Matrix([[self.a, self.b, 1]]) @ self.B

    @property
    def x(self):  # noqa: E743
        return Matrix([[self.c, self.d, 0]]) @ self.B

    @property
    def k(self):  # noqa: E743
        return Matrix([[self.e, self.f, 0]]) @ self.B

    @property
    def y(self):  # noqa: E743
        return Matrix([[0, 0, 1]]) @ self.B

    def constraint(self):
        return (self.x, self.k, self.y)

    def linicrypt_is_secure(self):
        a, b, c, d, e, f = self.a, self.b, self.c, self.d, self.e, self.f
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


def H2_permute_constraints():
    # Iterate over all combinations of a, b, c, d, e, f in {0,1}
    print("Permutation attack c1 <-> c2:")
    print("-----------------------------")
    count = 0
    for a, b, c, d, e, f in product([0, 1], repeat=6):
        f_compression = PGVComporessionFunction(a, b, c, d, e, f)
        H_f = MerkleDamgard(f_compression, 2)

        B = Matrix(H_f.dimension, H_f.dimension, symbols(f"B0:{H_f.dimension ** 2}"))
        permute_eqs = H_f.permute_eqs(1, 2, B)
        out_inv = H_f.output_invariant(B)
        equations = permute_eqs + [out_inv]
        solution = solve(equations, B)

        if solution:
            B_substituted = B.subs(solution)
            if B_substituted != eye(5):
                count += 1
                print(f_compression)
                pprint(B_substituted)
                print()
        else:
            pass
    print(count)


def H2_collapse_constraints(force_collision=True):
    # Iterate over all combinations of a, b, c, d, e, f in {0,1}
    print("Collapse attack c1 <-| c2:")
    print("--------------------------")
    count = 0
    for a, b, c, d, e, f in product([0, 1], repeat=6):
        f_compression = PGVComporessionFunction(a, b, c, d, e, f)
        H_f = MerkleDamgard(f_compression, 2)

        B = Matrix(H_f.dimension, H_f.dimension, symbols(f"B0:{H_f.dimension ** 2}"))
        equations = H_f.collapse_eqs(1, 2, B)
        if force_collision:
            equations += [H_f.output_h0(B)]
        solution = solve(equations, B)

        if solution:
            B_substituted = B.subs(solution)
            if B_substituted != eye(5):
                count += 1
                print(f_compression)
                pprint(B_substituted)
                if (
                    f_compression.brs_category() == "g"
                    and f_compression.pgv_category()[0] == "B"
                ):
                    pprint(B_substituted.eigenvects())
                print()
        else:
            pass
    print(count)


def Hn_cycle_constraints(n=3, force_collision=True):
    # Iterate over all combinations of a, b, c, d, e, f in {0,1}
    print("Permutation attack c1 -> c2 -> ... -> cn -> c1:")
    print("-----------------------------------------------")
    count = 0
    for a, b, c, d, e, f in product([0, 1], repeat=6):
        f_compression = PGVComporessionFunction(a, b, c, d, e, f)
        H_f = MerkleDamgard(f_compression, n)

        B = Matrix(H_f.dimension, H_f.dimension, symbols(f"B0:{H_f.dimension ** 2}"))
        equations = H_f.cycle_eqs(B)
        if force_collision:
            out_inv = H_f.output_invariant(B)
            equations += [out_inv]
        else:
            iv_inv = H_f.iv_invariant(B)
            equations += [iv_inv]
        solution = solve(equations, B)

        if solution:
            B_substituted = B.subs(solution)
            if B_substituted != eye(H_f.dimension):
                count += 1
                print(f_compression)
                pprint(B_substituted)
                if (
                    f_compression.brs_category() == "g"
                    and f_compression.pgv_category()[0] == "B"
                ):
                    pprint(B_substituted.eigenvects())
                    pprint(B_substituted.transpose().eigenvects())
                print()
        else:
            pass
    print(count)


# This doesnt make sense yet, but there has to be something there...
# def H2_flip_permute_constraints():
#     # Iterate over all combinations of a, b, c, d, e, f in {0,1}
#     print("Permutation attack on c_i = x <-k-> y to y <-k-> x:")
#     for a, b, c, d, e, f in product([0, 1], repeat=6):
#         f_compression = PGVComporessionFunction(a, b, c, d, e, f)
#         H_f = MerkleDamgard(f_compression, 2)

#         B = Matrix(H_f.dimension, H_f.dimension, symbols(f"B0:{H_f.dimension ** 2}"))
#         flip_1 = H_f.flip_permute(1, B)
#         out_inv = H_f.ouptut_invariant(B)
#         equations = flip_1 + [out_inv]
#         solution = solve(equations, B)

#         if solution:
#             print(f_compression)
#             B_substituted = B.subs(solution)
#             pprint(B_substituted)
#         else:
#             pass

if __name__ == "__main__":
    # H2_permute_constraints()
    Hn_cycle_constraints(n=4, force_collision=False)
    # H2_collapse_constraints(force_collision=True)
    # H2_flip_permute_constraints()
