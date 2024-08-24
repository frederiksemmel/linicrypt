from itertools import pairwise, permutations
from typing import Iterator

import numpy as np
from galois import FieldArray
from loguru import logger
from more_itertools import set_partitions
from tqdm import tqdm

from linicrypt_solver import Constraint, DualVector
from linicrypt_solver.field import GF
from linicrypt_solver.ideal_cipher import ConstraintE
from linicrypt_solver.random_oracle import ConstraintH
from linicrypt_solver.utils import stack_matrices


Partition = list[list[int]]


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

    def is_solution_ordering(self, fixing: FieldArray) -> bool:
        if len(self.cs) == 0:
            return True
        dim = self.dim()
        assert fixing.shape[1] == dim
        for i, c in enumerate(self.cs):
            if not c.is_solvable(fixing):
                logger.debug(f"Not solution ordering because of {i}: {c}")
                return False
            fixing = stack_matrices(fixing, c.fixing_matrix()).row_space()
        return True

    def is_solvable_brute_force(self, fixing: FieldArray) -> bool:
        for permuted_cs in permutations(self.cs):
            C_permuted = Constraints(list(permuted_cs))
            if C_permuted.is_solution_ordering(fixing):
                logger.info(f"Found solution ordering {permuted_cs}")
                return True
        return False

    def is_solvable(self, fixing: FieldArray) -> bool:
        ordering = []
        remaining = self.cs  # we are not modifying remaing, so this is ok
        # We go through self.cs and choose a constraint that is solvable
        # Then we repeat the process until we have an ordering
        # If an ordering exists, then this process will find an ordering (might be a different one)
        # Need to prove this across different constraint types

        # This while loop will finish in <= n loops
        while len(remaining) > 0:
            # TODO Use a function from itertools or more_itertools to partition remaining into a singleton set and the rest
            for i in range(len(remaining)):
                c = remaining[i]
                rest = remaining[:i] + remaining[i + 1 :]
                # Check if the constraint can be solved given the fixing of the rest
                # If it can, add it to the ordering and update the remaining constraints
                # If it cannot, continue the loop
                fixing_rest = GF(
                    np.concatenate([c.fixing_matrix() for c in rest] + [fixing])
                )
                logger.debug(f"c=\n{c}")
                logger.debug(f"rest=\n{rest}")
                logger.debug(f"fixing_rest=\n{fixing_rest}")
                if c.is_solvable(fixing_rest):
                    logger.debug(f"solving remaining {len(remaining)}: {c} is solvable")
                    ordering = [c] + ordering
                    remaining = rest
                    break
            # here we have found no solvable constraint, so the whole set has to be unsolvable
            else:
                logger.debug(f"solving remaining {len(remaining)}: nothing is solvable")
                return False

        # If we completed the while loop, ordering is a solution ordering
        assert Constraints(ordering).is_solution_ordering(fixing)
        return True

    def find_solvable_subspaces(
        self, fixing: FieldArray | None = None
    ) -> Iterator[tuple[Partition, FieldArray]]:
        if fixing is None:
            fixing = GF.Zeros((1, self.dim()))
        n = len(self.cs)

        # subspaces = []

        # https://codegolf.stackexchange.com/questions/132379/output-the-n-th-bell-number
        # https://en.wikipedia.org/wiki/Partition_of_a_set
        def bell_number(n, k=0):
            return n < 1 or k * bell_number(n - 1, k) + bell_number(n - 1, k + 1)

        # todo len
        for partition in tqdm(set_partitions(range(n)), total=bell_number(n)):
            logger.info(f"collapsing {partition}")
            collapsed_C, subspace = self.collapse(partition)
            collapsed_fixing = fixing @ subspace
            if collapsed_C.is_proper() and collapsed_C.is_solvable(collapsed_fixing):
                yield (partition, subspace)

    def find_solvable_subspaces_outside(
        self, W: FieldArray, fixing: FieldArray | None = None
    ) -> Iterator[tuple[Partition, FieldArray]]:
        dim_W = len(W.column_space())

        def is_outside_W(subspace):
            W_plus = stack_matrices(W, subspace, axis=1).column_space()
            assert len(W_plus) >= dim_W
            return len(W_plus) > dim_W

        return filter(
            lambda t: is_outside_W(t[1]), self.find_solvable_subspaces(fixing)
        )

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
                diff = stack_matrices(diff, self.cs[i].difference_matrix(self.cs[j]))

        logger.debug(f"diff matrix:\n{diff}")

        f_matrix = diff.null_space().transpose()
        logger.debug(f"collapsing space:\n{f_matrix}")
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
