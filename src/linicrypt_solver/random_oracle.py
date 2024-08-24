from typing import Self
import numpy as np
from galois import FieldArray
from loguru import logger

from linicrypt_solver.field import GF
from linicrypt_solver import Constraint, DualVector
from linicrypt_solver.utils import stack_matrices


class ConstraintH(Constraint):
    def __init__(self, q: DualVector, a: DualVector):
        if isinstance(q, list):
            q = np.array([q])
        if isinstance(a, list):
            a = np.array([a])
        m, n = a.shape
        assert m == 1
        assert n == q.shape[1]
        self.q = GF(q)
        self.a = GF(a)

    def fixing_matrix(self) -> FieldArray:
        return stack_matrices(self.q, self.a)

    def map(self, f: FieldArray) -> "ConstraintH":
        q = self.q @ f
        a = self.a @ f
        return ConstraintH(q, a)

    def dim(self):
        assert self.q.shape[1] == self.a.shape[1]
        return self.a.shape[1]

    # Returns the new fixed space
    def is_solvable(self, fixing: FieldArray) -> bool:
        fixing = stack_matrices(fixing, self.q).row_space()
        new_fixing_space = stack_matrices(fixing, self.a).row_space()
        assert len(fixing) <= len(new_fixing_space)
        if len(fixing.row_space()) == len(new_fixing_space.row_space()):
            logger.debug(f"{self.a} is contained in:\n{fixing}")
            return False
        return True

    def is_proper(self, fixed_constraints: list[Self]) -> bool:
        return all((c.q != self.q).any() for c in fixed_constraints)

    def __repr__(self):
        return f"{self.q[0]} |-> {self.a[0]}"

    def __eq__(self, other) -> bool:
        return (self.q == other.q).all() and (self.a == other.a).all()
