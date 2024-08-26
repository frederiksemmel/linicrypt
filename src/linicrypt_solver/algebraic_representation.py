import numpy as np
from galois import FieldArray
from loguru import logger
from typing import Iterator

from enum import Enum

from linicrypt_solver.field import GF
from linicrypt_solver.solvable import Constraints, Partition
from linicrypt_solver.utils import stack_matrices, embed_left, embed_right

SimpleAttack = tuple[Partition, FieldArray, Constraints]


class PartitionCompare(Enum):
    FINER = "finer"
    COARSER = "coarser"
    EQUAL = "equal"
    UNCOMPARABLE = "uncomparable"


def compare_partitions(partition_a, partition_b) -> PartitionCompare:
    a_finer_b = is_finer(partition_a, partition_b)
    b_finer_a = is_finer(partition_b, partition_a)
    if a_finer_b and b_finer_a:
        return PartitionCompare.EQUAL
    elif a_finer_b:
        return PartitionCompare.FINER
    elif b_finer_a:
        return PartitionCompare.COARSER
    else:
        return PartitionCompare.UNCOMPARABLE


def is_subset(set_a, set_b):
    return all(x in set_b for x in set_a)


def is_finer(partition_a, partition_b) -> bool:
    return all(
        any(is_subset(set_a, set_b) for set_b in partition_b) for set_a in partition_a
    )


def test_subset():
    a = [2, 3]
    b = [1, 2, 3]
    assert is_subset(a, b)


def test_coarser():
    a = [[1, 2], [3]]
    b = [[1], [2], [3]]
    assert compare_partitions(a, b) == PartitionCompare.COARSER


def test_finer():
    a = [[1, 2], [3]]
    b = [[1], [2], [3]]
    assert compare_partitions(b, a) == PartitionCompare.FINER


def test_uncomparable():
    a = [[1, 2], [3]]
    b = [[1], [2, 3]]
    assert compare_partitions(a, b) == PartitionCompare.UNCOMPARABLE
    assert compare_partitions(b, a) == PartitionCompare.UNCOMPARABLE


def test_equal():
    a = [[1, 2], [3]]
    b = [[1, 2], [3]]
    assert compare_partitions(a, b) == PartitionCompare.EQUAL
    assert compare_partitions(b, a) == PartitionCompare.EQUAL


class Attack:
    def __init__(
        self,
        partition: Partition,
        subspace: FieldArray | None,
        fixing: FieldArray | None,
        solution: Constraints,
    ):
        self.original_partition = partition
        self.frozen_partition = frozenset(frozenset(subset) for subset in partition)
        self.subspace = subspace
        self.fixing = fixing
        self.solution = solution

    @property
    def partition(self):
        return self.original_partition

    def __repr__(self) -> str:
        lines = [
            f"partition:\n{self.original_partition}",
            f"fixing:\n{self.fixing}",
            f"subspace:\n{self.subspace}",
            f"solution:\n{self.solution}",
        ]
        return "\n".join(lines)

    def __key(self):
        return self.frozen_partition

    def __hash__(self):
        return hash(self.__key())

    def __eq__(self, other):
        if isinstance(other, Attack):
            return self.__key() == other.__key()
        return NotImplemented


def maximal_attacks(attacks: Iterator[Attack]) -> list[Attack]:
    # If the partition is finer, the subspace is larger
    # partitions are only partially ordered, so we might have multiple maxima
    maxima: list[Attack] = []
    for attack in attacks:
        part = attack.partition
        part_is_uncomparable = True
        new_maxima = set()
        for maximum_attack in maxima:
            maximum_part = maximum_attack.partition
            comparision = compare_partitions(part, maximum_part)
            if comparision == PartitionCompare.FINER:
                logger.debug(f"{part} is finer than some in {maxima}")
                part_is_uncomparable = False
                new_maxima.add(attack)
            elif comparision == PartitionCompare.COARSER:
                logger.debug(f"{part} is coarser than some in {maxima}")
                part_is_uncomparable = False
                new_maxima.add(maximum_attack)
            elif comparision == PartitionCompare.EQUAL:
                logger.debug(f"{part} is equal to some in {maxima}")
                part_is_uncomparable = False
                new_maxima.add(maximum_attack)
            elif comparision == PartitionCompare.UNCOMPARABLE:
                logger.debug(f"{part} is uncomparable to some in {maxima}")
                new_maxima.add(maximum_attack)

        if part_is_uncomparable:
            logger.debug(f"{part} is umcomparable to all {maxima}")
            new_maxima.add(attack)

        maxima = list(new_maxima)
    return list(maxima)


def test_maximal_attack():
    zeros = GF.Zeros((1, 1))
    cs = Constraints([])
    a = Attack([[1, 2], [3]], zeros, zeros, cs)
    b = Attack([[1], [2, 3]], zeros, zeros, cs)
    assert maximal_attacks(iter([a, b])) == [a, b]

    a = Attack([[1, 2], [3]], zeros, zeros, cs)
    b = Attack([[1], [2], [3]], zeros, zeros, cs)
    assert maximal_attacks(iter([a, b])) == [b]

    a = Attack([[1, 2], [3]], zeros, zeros, cs)
    b = Attack([[1], [2], [3]], zeros, zeros, cs)
    c = Attack([[1], [2, 3]], zeros, zeros, cs)
    assert maximal_attacks(iter([a, b, c])) == [b]

    a = Attack([[1, 2], [3]], zeros, zeros, cs)
    b = Attack([[1], [2, 3]], zeros, zeros, cs)
    c = Attack([[1, 2, 3]], zeros, zeros, cs)
    assert maximal_attacks(iter([a, b, c])) == [a, b]

    a = Attack([[2], [1, 3], [0]], zeros, zeros, cs)
    b = Attack([[2], [1, 3], [0]], zeros, zeros, cs)
    assert maximal_attacks(iter([a, b])) == [a]


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

    def collapse_output_f(self):
        dim = self.dim()
        S = stack_matrices(GF.Identity(dim), GF.Identity(dim))
        ker_output = self.output.null_space().transpose()
        left_ker_output = stack_matrices(
            ker_output, GF.Zeros((dim, ker_output.shape[1]))
        )
        f = stack_matrices(S, left_ker_output, axis=1)
        return f

    def all_maximal_collision_attacks(self) -> list[Attack]:
        return maximal_attacks(self.list_collision_attacks())

    def list_collision_attacks(self) -> Iterator[Attack]:
        S = stack_matrices(GF.Identity(self.dim()), GF.Identity(self.dim()))
        C_join = self.cs.construct_joined()
        dim = C_join.dim()
        output_collapse = embed_left(self.output, dim) - embed_right(self.output, dim)
        # f = output_collapse.null_space().transpose()
        f = self.collapse_output_f()
        # print(f)
        # print(output_collapse)
        assert (output_collapse @ f == GF.Zeros((1, 1))).all()

        # Robust way to compute the preimage of S
        # Annihlator of S called S^0 are the dual vectors that are zero on S
        S_0 = S.left_null_space()
        # This is f^*(S^0)
        f_pullback_S_0 = S_0 @ f
        # We have f^-1(S) = f^*(S^0))^0. Because S is in the image of f, this f(f^-1(S)) = S
        preimage_S = f_pullback_S_0.null_space().transpose()
        # This could be wrong, because preimage_S might be a different basis
        # of the same subspace. But it seems the left_null_space and null_space
        # algorithms of the galois package are such that this actually works.
        # So each column of preimage_S is actually the preimage of each column of S
        assert (f @ preimage_S == S).all()

        def is_outside_S(subspace):
            W_plus = stack_matrices(S, subspace, axis=1).column_space()
            assert len(W_plus) >= self.dim()
            return len(W_plus) > self.dim()

        C_joined_f = C_join.map(f)
        subspaces_iter = C_joined_f.find_solvable_subspaces_outside(preimage_S)
        for part, subspace in subspaces_iter:
            logger.info("Found solvable subspace:")
            logger.info(f"Partition of the constraints is {part}")
            logger.info("Solvable subspace of F^(2d) is")
            logger.info(f @ subspace)
            logger.info("Solvable constraints in that subspace are")
            collapsed_constraints = C_join.map(f @ subspace)
            zero = GF.Zeros((1, collapsed_constraints.dim()))
            solution = collapsed_constraints.find_solution_ordering(fixing=zero)
            logger.info(solution)
            assert solution is not None
            assert is_outside_S(f @ subspace)
            yield Attack(part, f @ subspace, None, C_join.map(f @ subspace))

    def is_collision_resistant(self):
        return any(True for _ in self.list_collision_attacks())

    def list_second_preimage_attacks(self):
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
            yield (part, f @ subspace, C_join.map(f @ subspace))

    def is_second_preimage_resistant(self):
        return any(True for _ in self.list_second_preimage_attacks())
