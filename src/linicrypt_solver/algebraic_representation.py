import numpy as np
from galois import FieldArray
from loguru import logger

from linicrypt_solver.field import GF
from linicrypt_solver.solvable import Constraints
from linicrypt_solver.utils import stack_matrices, embed_left, embed_right


class AlgebraicRep:
    def __init__(self, cs: Constraints, fixing: FieldArray, output: FieldArray):
        dim_cs = cs.dim()
        dim_fixing = fixing.shape[1]
        dim_output = output.shape[1]
        assert dim_cs == dim_fixing and dim_fixing == dim_output

        self.cs = cs
        self.fixing = fixing
        self.output = output

    def map(self, f: FieldArray) -> "AlgebraicRep":
        fixing = (self.fixing @ f).row_space()
        output = self.output @ f
        return AlgebraicRep(self.cs.map(f), fixing, output)

    def dim(self) -> int:
        return self.cs.dim()

    def merge(self, other: "AlgebraicRep"):
        assert self.dim() == other.dim()
        self.fixing = stack_matrices(self.fixing, other.fixing).row_space()
        self.output = stack_matrices(self.output, other.output)
        self.cs.merge(other.cs)

    def embed_left(self, dim: int) -> "AlgebraicRep":
        own_dim = self.dim()
        assert dim >= own_dim
        identity = GF.Identity(own_dim)
        zeros = GF.Zeros((own_dim, dim - own_dim))
        f = GF(np.concatenate((identity, zeros), axis=1))
        return self.map(f)

    def embed_right(self, dim: int) -> "AlgebraicRep":
        own_dim = self.dim()
        assert dim >= own_dim
        identity = GF.Identity(own_dim)
        zeros = GF.Zeros((own_dim, dim - own_dim))
        f = GF(np.concatenate((zeros, identity), axis=1))
        return self.map(f)

    def __repr__(self) -> str:
        lines = [
            f"I:\n{self.fixing}",
            f"O:\n{self.output}",
            f"C:\n{self.cs}",
        ]
        return "\n".join(lines)

    def is_collision_resistant(self):
        S = stack_matrices(GF.Identity(self.dim()), GF.Identity(self.dim()))
        C_join = self.cs.construct_joined()
        dim = C_join.dim()
        output_collapse = embed_left(self.output, dim) - embed_right(self.output, dim)
        f = output_collapse.null_space().transpose()

        # Robust way to compute the preimage of S
        # Annihlator of S called S^0 are the dual vectors that are zero on S
        S_0 = S.left_null_space()
        # This is f^*(S^0)
        f_pullback_S_0 = S_0 @ f
        # We have f^-1(S) = f^*(S^0))^0. Because S is in the image of f, this f(f^-1(S)) = S
        preimage_S = f_pullback_S_0.null_space().transpose()
        assert (f @ preimage_S == S).all()

        def is_outside_S(subspace):
            W_plus = stack_matrices(S, subspace, axis=1).column_space()
            assert len(W_plus) >= self.dim()
            return len(W_plus) > self.dim()

        C_joined_f = C_join.map(f)
        subspaces_iter = C_joined_f.find_solvable_subspaces_outside(preimage_S)
        for part, subspace in subspaces_iter:
            logger.info("Found solvable subspace:")
            logger.info(part)
            logger.info(subspace)
            solution = C_joined_f.map(subspace)
            logger.info(solution)
            dim = solution.dim()
            logger.info(f @ subspace)
            logger.info(C_join.map(f @ subspace))
            assert is_outside_S(f @ subspace)
            # logger.info(C_joined_f.map(subspace).is_solvable(fixing=GF.Zeros((1, dim))))
            return False
        return True

    def is_second_preimage_resistant(self):
        S = stack_matrices(GF.Identity(self.dim()), GF.Identity(self.dim()))
        C_join = self.cs.construct_joined()
        dim = C_join.dim()
        output_collapse = embed_left(self.output, dim) - embed_right(self.output, dim)
        f = output_collapse.null_space().transpose()

        # Robust way to compute the preimage of S
        # Annihlator of S called S^0 are the dual vectors that are zero on S
        S_0 = S.left_null_space()
        # This is f^*(S^0)
        f_pullback_S_0 = S_0 @ f
        # We have f^-1(S) = f^*(S^0))^0. Because S is in the image of f, this f(f^-1(S)) = S
        preimage_S = f_pullback_S_0.null_space().transpose()
        assert (f @ preimage_S == S).all()

        I_1 = embed_left(self.fixing, dim)
        logger.debug(f"Left input is:\n{I_1}")
        I_1_f = I_1 @ f
        logger.debug(f"Left input after pullback is:\n{I_1_f}")

        C_joined_f = C_join.map(f)
        subspaces_iter = C_joined_f.find_solvable_subspaces_outside(
            preimage_S, fixing=I_1_f
        )
        for part, subspace in subspaces_iter:
            logger.info(part)
            logger.info(subspace)
            return False
        return True
