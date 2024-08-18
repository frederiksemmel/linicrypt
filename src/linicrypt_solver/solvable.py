import sys
from itertools import pairwise, permutations

import galois
import numpy as np
from galois import FieldArray
from loguru import logger
from more_itertools import set_partitions

# Configure logger to print the log message on a new line
logger.remove()  # Remove the default handler
logger.add(
    sink=sys.stderr,
    level="INFO",
    format="{time:YYYY-MM-DD HH:mm:ss} | {level}\n{message}",  # Add newline before {message}
)

field_size = 2**4
GF = galois.GF(field_size)


def stack_matrices(A: FieldArray, B: FieldArray, axis=0) -> FieldArray:
    return GF(np.concatenate((A, B), axis=axis))


def embed_left(A: FieldArray) -> FieldArray:
    zeros = GF.Zeros(A.shape)
    return GF(np.concatenate((A, zeros), axis=1))


def embed_right(A: FieldArray) -> FieldArray:
    zeros = GF.Zeros(A.shape)
    return GF(np.concatenate((zeros, A), axis=1))


class Constraint:
    def __init__(self, q: np.ndarray, a: np.ndarray):
        m, n = a.shape
        assert m == 1
        assert n == q.shape[1]

        self.q = GF(q)
        self.a = GF(a)

    def difference_matrix(self, other: "Constraint") -> FieldArray:
        q_row = self.q - other.q
        a_row = self.a - other.a
        return stack_matrices(q_row, a_row)

    def map(self, f: FieldArray) -> "Constraint":
        q = self.q @ f
        a = self.a @ f
        return Constraint(q, a)

    def dim(self):
        assert self.q.shape[1] == self.a.shape[1]
        return self.a.shape[1]

    def __repr__(self):
        return f"{self.q[0]} |-> {self.a[0]}"

    def __eq__(self, other) -> bool:
        return (self.q == other.q).all() and (self.a == other.a).all()

    def embed_left(self) -> "Constraint":
        q = embed_left(self.q)
        a = embed_left(self.a)
        return Constraint(q, a)

    def embed_right(self) -> "Constraint":
        q = embed_right(self.q)
        a = embed_right(self.a)
        return Constraint(q, a)


class Constraints:
    def __init__(self, cs: list[Constraint]):
        ordered_set = []
        for c in cs:
            if c not in ordered_set:
                ordered_set.append(c)
        self.cs: list[Constraint] = ordered_set

    @staticmethod
    def from_repr(representation: list[tuple[list[int], list[int]]]) -> "Constraints":
        cs = []
        for q, a in representation:
            c = Constraint(np.array([q]), np.array([a]))
            cs.append(c)
        return Constraints(cs)

    def is_solution_ordering(self) -> bool:
        if len(self.cs) == 0:
            return True
        fixing_space = self.cs[0].q
        qs = []
        logger.debug(f"Checking solution ordering of {self.cs}")
        for i, c in enumerate(self.cs):
            fixing_space = stack_matrices(fixing_space, c.q).row_space()
            new_fixing_space = stack_matrices(fixing_space, c.a).row_space()
            assert len(fixing_space) <= len(new_fixing_space)
            if len(fixing_space.row_space()) == len(new_fixing_space.row_space()):
                logger.debug(f"Not solution ordering because of {i}: {c.a}")
                return False
            if any((c.q == q).all() for q in qs):
                logger.debug(
                    f"Not solution ordering because {c.q} repeated with {c.a} different"
                )
                return False
            qs.append(c.q)
            fixing_space = new_fixing_space

        logger.info(f"{self.cs} is solution ordering")
        return True

    def is_solvable(self) -> bool:
        for permuted_cs in permutations(self.cs):
            C_permuted = Constraints(list(permuted_cs))
            if C_permuted.is_solution_ordering():
                return True

        return False

    def find_solvable_subspaces(self) -> list[tuple[tuple[int, ...], FieldArray]]:
        n = len(self.cs)

        subspaces = []
        # todo len
        for partition in set_partitions(range(n)):
            logger.info(f"collapsing {partition}")
            collapsed_C, subspace = self.collapse(partition)
            if collapsed_C.is_solvable():
                subspaces.append((partition, subspace))
        return subspaces

    def find_solvable_subspaces_outside(
        self, W: FieldArray
    ) -> list[tuple[tuple[int, ...], FieldArray]]:
        dim_W = len(W.column_space())

        def is_outside_W(subspace):
            W_plus = stack_matrices(W, subspace, axis=1).column_space()
            assert len(W_plus) >= dim_W
            return len(W_plus) > dim_W

        return [
            (part, subspace)
            for part, subspace in self.find_solvable_subspaces()
            if is_outside_W(subspace)
        ]

    def collapse_pair(self, i: int, j: int) -> "Constraints":
        assert i != j
        assert i < len(self.cs)
        assert j < len(self.cs)

        diff = self.cs[i].difference_matrix(self.cs[j])
        f_matrix = diff.null_space().transpose()
        return self.map(f_matrix)

    def collapse(self, partition: list[list[int]]) -> tuple["Constraints", FieldArray]:
        # todo all different todo dont collapse two of the same todo
        assert all(i < len(self.cs) for collapse in partition for i in collapse)

        # todo what if no constraints
        d = self.cs[0].dim()
        diff = GF.Zeros((1, d))
        for collapse in partition:
            for i, j in pairwise(collapse):
                logger.debug(f"in collapse collapsing {i},{j}")
                diff = stack_matrices(diff, self.cs[i].difference_matrix(self.cs[j]))

        logger.debug(f"diff {diff}")

        f_matrix = diff.null_space().transpose()
        return (self.map(f_matrix), f_matrix)

    def map(self, f: FieldArray) -> "Constraints":
        return Constraints([c.map(f) for c in self.cs])

    def __repr__(self):
        lines = []
        for c in self.cs:
            lines.append((f"{c}"))
        return "\n".join(lines)

    def construct_joined(self):
        cs_1 = [c.embed_left() for c in self.cs]
        cs_2 = [c.embed_right() for c in self.cs]

        return Constraints(cs_1 + cs_2)


if __name__ == "__main__":
    C = Constraints.from_repr(
        [
            ([1, 0, 0, 0, 0], [0, 0, 1, 0, 0]),
            ([0, 0, 1, 0, 0], [0, 0, 0, 1, 0]),
            ([0, 1, 0, 0, 0], [0, 0, 0, 0, 1]),
        ]
    )
    output = GF([[1, 1, 1, 1, 1]])
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
    output_collapse = embed_left(output) - embed_right(output)
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
