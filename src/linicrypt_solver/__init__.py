from abc import ABC, abstractmethod
from typing import Self

import numpy as np
from galois import FieldArray

from linicrypt_solver.field import GF


# TODO flesh out interface and add make constraints work with both types of constraints simultaneously
class Constraint(ABC):
    @abstractmethod
    def difference_matrix(self, other: Self) -> FieldArray:
        pass

    @abstractmethod
    def dim(self) -> int:
        pass

    @abstractmethod
    def map(self, f: FieldArray) -> Self:
        pass

    @abstractmethod
    def is_solvable(self, fixing: FieldArray) -> None | FieldArray:
        pass

    @abstractmethod
    def is_proper(self, fixed_constraints: list[Self]) -> bool:
        pass

    @abstractmethod
    def __eq__(self, other) -> bool:
        pass

    @abstractmethod
    def __repr__(self) -> str:
        pass

    def embed_left(self, dim: int) -> Self:
        own_dim = self.dim()
        assert dim >= own_dim
        identity = GF.Identity(own_dim)
        zeros = GF.Zeros((own_dim, dim - own_dim))
        f = GF(np.concatenate((identity, zeros), axis=1))
        return self.map(f)

    def embed_right(self, dim: int) -> Self:
        own_dim = self.dim()
        assert dim >= own_dim
        identity = GF.Identity(own_dim)
        zeros = GF.Zeros((own_dim, dim - own_dim))
        f = GF(np.concatenate((zeros, identity), axis=1))
        return self.map(f)


DualVector = list[int] | FieldArray | np.ndarray
