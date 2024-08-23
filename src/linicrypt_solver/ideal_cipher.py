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

    def difference_matrix(self, other: "ConstraintE") -> FieldArray:
        x_diff = self.x - other.x
        k_diff = self.k - other.k
        y_diff = self.y - other.y
        return GF(np.concatenate((x_diff, k_diff, y_diff)))

    def map(self, f: FieldArray) -> "ConstraintE":
        x = self.x @ f
        k = self.k @ f
        y = self.y @ f
        return ConstraintE(x, k, y)

    def is_solvable_enc(self, fixing: FieldArray) -> None | FieldArray:
        fixing = GF(np.concatenate((fixing, self.x, self.k))).row_space()
        new_fixing_space = stack_matrices(fixing, self.y).row_space()
        assert len(fixing) <= len(new_fixing_space)
        if len(fixing.row_space()) == len(new_fixing_space.row_space()):
            logger.debug(
                f"solvable_enc: y = {self.y} is contained in:\n{fixing} + <x,y>"
            )
            return None
        return new_fixing_space

    def is_solvable_dec(self, fixing: FieldArray) -> None | FieldArray:
        fixing = GF(np.concatenate((fixing, self.k, self.y))).row_space()
        new_fixing_space = stack_matrices(fixing, self.x).row_space()
        assert len(fixing) <= len(new_fixing_space)
        if len(fixing.row_space()) == len(new_fixing_space.row_space()):
            logger.debug(
                f"solvable_dec: y = {self.y} is contained in:\n{fixing} + <x,y>"
            )
            return None
        return new_fixing_space

    def is_solvable(self, fixing: FieldArray) -> None | FieldArray:
        enc_space = self.is_solvable_enc(fixing)
        dec_space = self.is_solvable_dec(fixing)
        if enc_space is None and dec_space is None:
            return None
        return enc_space

    def is_proper(self, fixed_constraints: list[Self]) -> bool:
        xk_different = all((c.x, c.k) != (self.x, self.k) for c in fixed_constraints)
        ky_different = all((c.k, c.y) != (self.k, self.y) for c in fixed_constraints)
        return xk_different and ky_different

    def dim(self):
        return self.k.shape[1]

    def components(self) -> FieldArray:
        return GF(np.concatenate((self.x, self.x, self.y)))

    def __repr__(self):
        return f"{self.x[0]} <- {self.k[0]} -> {self.y[0]}"

    def __eq__(self, other) -> bool:
        assert isinstance(other, ConstraintE)
        return (self.components() == other.components()).all()
