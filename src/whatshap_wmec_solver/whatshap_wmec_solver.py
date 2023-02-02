#!/usr/bin/env python
# coding: utf-8

from __future__ import annotations
from typing import Iterable, Sequence
import numpy as np
import whatshap.core as wh
from whatshap.testhelpers import canonic_index_to_biallelic_gt


class wMECSolver:
    def __init__(self, readset):
        self.readset = readset

    @classmethod
    def from_matrix(cls, matrix, weights=None):
        rs = wh.ReadSet()
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
            rs.add(read)
        return cls(rs)

    @classmethod
    def from_matrix_file(cls, matrix_file_handle):
        matrix = []
        for line in matrix_file_handle:
            row = [int(c) if c.isdigit() else -1 for c in line]
            matrix.append(row)
        return cls.from_matrix(matrix)

    @classmethod
    def from_compressed_matrix_file(cls, matrix_file):
        pass

    def compute_haplotypes(self, all_heterozygous=False):
        positions = self.readset.get_positions()
        recombcost = [1] * len(positions)
        pedigree = wh.Pedigree(wh.NumericSampleIds())
        genotype_likelihoods = [
            None if all_heterozygous else wh.PhredGenotypeLikelihoods([0, 0, 0])
        ] * len(positions)
        pedigree.add_individual(
            "individual0",
            [canonic_index_to_biallelic_gt(1) for i in range(len(positions))],
            genotype_likelihoods,
        )

        dp_table = wh.PedigreeDPTable(
            self.readset, recombcost, pedigree, distrust_genotypes=not all_heterozygous
        )

        superreads, transmission_vector = dp_table.get_super_reads()
        assert len(set(transmission_vector)) == 1

        cost = dp_table.get_optimal_cost()
        partition = dp_table.get_optimal_partitioning()
        hap1 = [v.allele for v in superreads[0][0]]
        hap2 = [v.allele for v in superreads[0][1]]

        # Represent ambiguous values with -1 instead of 3
        hap1 = [x if x != 3 else -1 for x in hap1]
        hap2 = [x if x != 3 else -1 for x in hap2]
        
        return (hap1, hap2), partition, cost, superreads


def solve_wMEC(
    allele_matrix: Iterable[Iterable[int]],
    weights: Iterable[Iterable[int]] | None = None,
) -> tuple[Sequence[int], Sequence[int]]:
    solver = wMECSolver.from_matrix(allele_matrix, weights=weights)
    (haplotype_1, haplotype_2), *_ = solver.compute_haplotypes()
    if haplotype_1 > haplotype_2:
        haplotype_1, haplotype_2 = haplotype_2, haplotype_1
    return haplotype_1, haplotype_2
