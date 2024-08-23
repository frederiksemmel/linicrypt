from itertools import pairwise, permutations

from linicrypt_solver.field import GF
from galois import FieldArray
from loguru import logger
from more_itertools import set_partitions

from linicrypt_solver import Constraint, DualVector
from linicrypt_solver.ideal_cipher import ConstraintE
from linicrypt_solver.random_oracle import ConstraintH
from linicrypt_solver.utils import stack_matrices


class Constraints:
    def __init__(self, cs: list[Constraint]):
        ordered_set = []
        for c in cs:
            if c not in ordered_set:
                ordered_set.append(c)
        self.cs: list[Constraint] = ordered_set

    def add(self, c: Constraint):
        if len(self.cs) > 0:
            assert c.dim() == self.cs[0].dim()
        self.cs.append(c)

    def merge(self, other: "Constraints"):
        assert self.dim() == other.dim()
        self.cs += other.cs

    @staticmethod
    def from_repr(representation: list[tuple[DualVector, ...]]) -> "Constraints":
        cs = []
        for c_repr in representation:
            if len(c_repr) == 2:
                q, a = c_repr
                c = ConstraintH(q, a)
            elif len(c_repr) == 3:
                x, k, y = c_repr
                c = ConstraintE(x, k, y)
            else:
                raise ValueError(f"Unknown constraint representation {c_repr}")
            cs.append(c)
        return Constraints(cs)

    def is_proper(self) -> bool:
        # check if some queries are exactly the same.
        # Then they should have the same answer vector (and the constraints should have
        # collapsed during construction, or in math term, thanks to describing constraints as a set)
        for i, c in enumerate(self.cs):
            if not c.is_proper(self.cs[:i]):
                logger.debug(
                    f"Not solution ordering because {c} is not proper with {self.cs[:i]}"
                )
                return False
        return True

    def is_solution_ordering(self, fixing: FieldArray | None = None) -> bool:
        if len(self.cs) == 0:
            return True
        dim = self.dim()
        if fixing is None:
            fixing = GF.Zeros((1, dim))
        for i, c in enumerate(self.cs):
            new_fixing = c.is_solvable(fixing)
            if new_fixing is None:
                logger.debug(f"Not solution ordering because of {i}: {c}")
                return False
            fixing = new_fixing
        return True

    def is_solvable(self, fixing: FieldArray | None = None) -> bool:
        for permuted_cs in permutations(self.cs):
            C_permuted = Constraints(list(permuted_cs))
            if C_permuted.is_solution_ordering(fixing):
                logger.info(f"Found solution ordering {permuted_cs}")
                return True
        return False

    def find_solvable_subspaces(
        self, fixing: FieldArray | None = None
    ) -> list[tuple[tuple[int, ...], FieldArray]]:
        n = len(self.cs)

        subspaces = []
        # todo len
        for partition in set_partitions(range(n)):
            logger.debug(f"collapsing {partition}")
            collapsed_C, subspace = self.collapse(partition)
            if collapsed_C.is_solvable(fixing):
                subspaces.append((partition, subspace))
        return subspaces

    def find_solvable_subspaces_outside(
        self, W: FieldArray, fixing: FieldArray | None = None
    ) -> list[tuple[tuple[int, ...], FieldArray]]:
        dim_W = len(W.column_space())

        def is_outside_W(subspace):
            W_plus = stack_matrices(W, subspace, axis=1).column_space()
            assert len(W_plus) >= dim_W
            return len(W_plus) > dim_W

        return [
            (part, subspace)
            for part, subspace in self.find_solvable_subspaces(fixing)
            if is_outside_W(subspace)
        ]

    def collapse_pair(self, i: int, j: int) -> "Constraints":
        assert i != j
        assert i < len(self.cs)
        assert j < len(self.cs)
        assert type(self.cs[i]) is type(self.cs[j])

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

        logger.debug(f"diff matrix:\n{diff}")

        f_matrix = diff.null_space().transpose()
        return (self.map(f_matrix), f_matrix)

    def map(self, f: FieldArray) -> "Constraints":
        return Constraints([c.map(f) for c in self.cs])

    def __repr__(self):
        lines = []
        for c in self.cs:
            lines.append((f"{c}"))
        return "\n".join(lines)

    def dim(self) -> int:
        if len(self.cs) == 0:
            return 0
        dim = self.cs[0].dim()
        assert all(dim == c.dim() for c in self.cs)
        return dim

    def construct_joined(self):
        dim = self.dim()
        cs_1 = [c.embed_left(dim * 2) for c in self.cs]
        cs_2 = [c.embed_right(dim * 2) for c in self.cs]

        return Constraints(cs_1 + cs_2)

    def embed_left(self, dim: int):
        return Constraints([c.embed_left(dim) for c in self.cs])

    def embed_right(self, dim: int):
        return Constraints([c.embed_right(dim) for c in self.cs])
