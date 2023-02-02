#!/usr/bin/env python
# coding: utf-8

from __future__ import annotations
from typing import Sequence
from dataclasses import dataclass
import whatshap.core as wh
from whatshap.testhelpers import canonic_index_to_biallelic_gt


# Type aliases
AlleleMatrix = Sequence[Sequence[int]]
WeightMatrix = Sequence[Sequence[float]]


@dataclass
class _wMECSolverResult:
    haplotypes: tuple[Sequence[int], Sequence[int]]
    partition: Sequence[int]
    cost: float


class wMECSolver:
    def __init__(self, matrix: AlleleMatrix, weights: WeightMatrix | None = None):
        # Sort reads and record indices
        sorted_indices = self._argsort_reads(matrix)
        _matrix = [matrix[i] for i in sorted_indices]
        _weights: WeightMatrix | None
        if weights is None:
            _weights = None
        else:
            _weights = [weights[i] for i in sorted_indices]
        self._index_mapping: dict[int, int] = {
            input_index: output_index
            for output_index, input_index in enumerate(sorted_indices)
        }

        readset = wh.ReadSet()
        mapping_quality = 50
        source_id = 0
        for i, row in enumerate(_matrix):
            assert len(row) > 1
            read_id = f"R{i+1}"
            read = wh.Read(read_id, mapping_quality, source_id)
            for j, column in enumerate(row):
                allele = column
                if allele == -1:
                    continue
                position = (j + 1) * 10  # Not sure why. See testhelpers.py:28
                w = _weights[i][j] if _weights is not None else 1
                read.add_variant(position, allele, w)
            readset.add(read)
        self._readset = readset

    @staticmethod
    def _argsort_reads(matrix: AlleleMatrix) -> Sequence[int]:
        # Reads need to be sorted by the position of the first non-ambiguous variant.
        # See https://github.com/whatshap/whatshap/blob/9882248c722b1020321fea6e9491c1cc5b75354b/src/columniterator.cpp (line #28)
        def get_first_variant_position(read: Sequence[int]) -> int:
            for i, allele in enumerate(read):
                if allele >= 0:
                    return i
            else:
                raise ValueError(f"Read contains only ambiguous alleles: {read}")

        read_indices = [(read, i) for i, read in enumerate(matrix)]
        sorted_reads = list(
            sorted(
                [(read, i) for i, read in enumerate(matrix)],
                key=lambda pair: get_first_variant_position(pair[0]),
            )
        )
        sorted_indices = [pair[1] for pair in sorted_reads]
        return sorted_indices

    def _reorder_partition(self, partition: Sequence[int]) -> Sequence[int]:
        reordered_partition: list[int] = []
        for input_index in sorted(self._index_mapping.keys()):
            output_index = self._index_mapping[input_index]
            reordered_partition.append(partition[output_index])
        return tuple(reordered_partition)

    def solve(self, allow_homozygousity: bool = True) -> _wMECSolverResult:
        positions = self._readset.get_positions()
        recombcost = [1] * len(positions)
        pedigree = wh.Pedigree(wh.NumericSampleIds())
        genotype_likelihoods = [
            None if not allow_homozygousity else wh.PhredGenotypeLikelihoods([0, 0, 0])
        ] * len(positions)
        pedigree.add_individual(
            "individual0",
            [canonic_index_to_biallelic_gt(1) for i in range(len(positions))],
            genotype_likelihoods,
        )

        dp_table = wh.PedigreeDPTable(
            self._readset, recombcost, pedigree, distrust_genotypes=allow_homozygousity
        )

        superreads, transmission_vector = dp_table.get_super_reads()
        assert len(set(transmission_vector)) == 1

        cost = dp_table.get_optimal_cost()
        partition = dp_table.get_optimal_partitioning()
        hap1: list[int] = [v.allele for v in superreads[0][0]]
        hap2: list[int] = [v.allele for v in superreads[0][1]]

        # Represent ambiguous values with -1 instead of 3
        hap1 = [x if x != 3 else -1 for x in hap1]
        hap2 = [x if x != 3 else -1 for x in hap2]

        return _wMECSolverResult(
            haplotypes=(tuple(hap1), tuple(hap2)),
            partition=self._reorder_partition(partition),
            cost=cost,
        )


def solve_wMEC(
    allele_matrix: AlleleMatrix,
    weights: WeightMatrix | None = None,
    *,
    allow_homozygousity=True,
) -> tuple[Sequence[int], Sequence[int]]:
    solver = wMECSolver(allele_matrix, weights=weights)
    result = solver.solve(allow_homozygousity=allow_homozygousity)
    haplotype_1, haplotype_2 = result.haplotypes
    return haplotype_1, haplotype_2
