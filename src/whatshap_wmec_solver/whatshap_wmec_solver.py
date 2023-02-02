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
        readset = wh.ReadSet()
        mapping_quality = 50
        source_id = 0
        for i, row in enumerate(matrix):
            assert len(row) > 1
            read_id = f"R{i+1}"
            read = wh.Read(read_id, mapping_quality, source_id)
            for j, column in enumerate(row):
                allele = column
                if allele == -1:
                    continue
                position = (j + 1) * 10  # Not sure why. See testhelpers.py:28
                w = weights[i][j] if weights else 1
                read.add_variant(position, allele, w)
            readset.add(read)
        self._readset = readset

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
            haplotypes=(tuple(hap1), tuple(hap2)), partition=tuple(partition), cost=cost
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
