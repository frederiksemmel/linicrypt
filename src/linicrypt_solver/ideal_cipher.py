import sys
from typing import Self

import numpy as np
from galois import FieldArray
from loguru import logger

from linicrypt_solver.field import GF
from linicrypt_solver.utils import stack_matrices
from linicrypt_solver import Constraint, DualVector

# Configure logger to print the log message on a new line
logger.remove()  # Remove the default handler
logger.add(
    sink=sys.stderr,
    level="INFO",
    format="{time:YYYY-MM-DD HH:mm:ss} | {level}\n{message}",  # Add newline before {message}
)


class ConstraintE(Constraint):
    def __init__(self, x_raw: DualVector, k_raw: DualVector, y_raw: DualVector):
        if isinstance(x_raw, list):
            x_raw = np.array([x_raw])
        if isinstance(k_raw, list):
            k_raw = np.array([k_raw])
        if isinstance(y_raw, list):
            y_raw = np.array([y_raw])
        x, k, y = GF(x_raw), GF(k_raw), GF(y_raw)
        assert x.shape == k.shape and k.shape == y.shape
        self.x = x
        self.k = k
        self.y = y

    def fixing_matrix(self) -> FieldArray:
        return GF(np.concatenate((self.x, self.k, self.y)))

    def map(self, f: FieldArray) -> "ConstraintE":
        x = self.x @ f
        k = self.k @ f
        y = self.y @ f
        return ConstraintE(x, k, y)

    def is_solvable_enc(self, fixing: FieldArray) -> bool:
        fixing = GF(np.concatenate((fixing, self.x, self.k))).row_space()
        new_fixing_space = stack_matrices(fixing, self.y).row_space()
        assert len(fixing) <= len(new_fixing_space)
        if len(fixing.row_space()) == len(new_fixing_space.row_space()):
            logger.debug(
                f"solvable_enc: y = {self.y} is contained in:\n{fixing} + <x,k>"
            )
            return False
        return True

    def is_solvable_dec(self, fixing: FieldArray) -> bool:
        fixing = GF(np.concatenate((fixing, self.k, self.y))).row_space()
        new_fixing_space = stack_matrices(fixing, self.x).row_space()
        assert len(fixing) <= len(new_fixing_space)
        if len(fixing.row_space()) == len(new_fixing_space.row_space()):
            logger.debug(
                f"solvable_dec: x = {self.x} is contained in:\n{fixing} + <k,y>"
            )
            return False
        return True

    def is_solvable(self, fixing: FieldArray) -> bool:
        return self.is_solvable_enc(fixing) or self.is_solvable_dec(fixing)

    def is_proper(self, fixed_constraints: list["ConstraintE"]) -> bool:
        my_matrix = self.fixing_matrix()
        my_xk = my_matrix[:2]
        my_ky = my_matrix[1:]
        for c in fixed_constraints:
            c_matrix = c.fixing_matrix()
            xk = c_matrix[:2]
            ky = c_matrix[1:]
            if (my_xk == xk).all() or (my_ky == ky).all():
                return False
        return True

    def dim(self):
        return self.k.shape[1]

    def components(self) -> FieldArray:
        return GF(np.concatenate((self.x, self.k, self.y)))

    def __repr__(self):
        return f"{self.x[0]} <- {self.k[0]} -> {self.y[0]}"
