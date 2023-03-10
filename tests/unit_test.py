#!/usr/bin/env python
# coding: utf-8

from whatshap_wmec_solver import wMECSolver, solve_wMEC


def test_simple_matrix():
    # Two reads with no conflicts
    allele_matrix = [[0, 0, 1, 1], [1, 1, 0, 0]]
    haplotypes = solve_wMEC(allele_matrix)
    assert set(haplotypes) == {(0, 0, 1, 1), (1, 1, 0, 0)}

    # The second variant is homozygous
    allele_matrix = [[0, 1, 0], [0, 1, 0], [0, 1, 1]]
    haplotypes = solve_wMEC(allele_matrix)
    assert set(haplotypes) == {(0, 1, 0), (0, 1, 1)}

    # The last read contains an error at the second variant
    allele_matrix = [[0, 1, 0], [1, 0, 1], [1, 0, 1], [1, 1, 1]]
    haplotypes = solve_wMEC(allele_matrix)
    assert set(haplotypes) == {(0, 1, 0), (1, 0, 1)}


def test_weights():
    # The last read conflicts with other reads at the second variant but has higher weights
    allele_matrix = [[0, 1, 0], [1, 0, 1], [1, 0, 1], [1, 1, 1]]
    weights = [[1, 1, 1], [1, 1, 1], [1, 1, 1], [10, 10, 10]]
    haplotypes = solve_wMEC(allele_matrix, weights)
    assert set(haplotypes) == {(0, 1, 0), (1, 1, 1)}


def test_ambiguous_matrix():
    # Read #3 conflicts with Read #2 at Variant #2. Without additional evidence, Variant #2 is ambiguous.
    allele_matrix = [[0, 1, 0], [1, 0, 1], [1, 1, 1]]
    haplotypes = solve_wMEC(allele_matrix)
    assert set(haplotypes) == {(0, 1, 0), (1, -1, 1)}


def test_reads_with_holes():
    allele_matrix = [[0, 1, 0], [1, 0, 1], [1, -1, 1]]
    haplotypes = solve_wMEC(allele_matrix)
    assert set(haplotypes) == {(0, 1, 0), (1, 0, 1)}


def test_multiple_phase_blocks():
    allele_matrix = [[0, 1, -1, -1], [1, 0, -1, -1], [-1, -1, 1, 0], [-1, -1, 0, 1]]
    haplotypes = solve_wMEC(allele_matrix)
    assert set(haplotypes) == {(0, 1, 1, 0), (1, 0, 0, 1)}


def test_allow_homozygous():
    allele_matrix = [[0, 1, 1], [0, 1, 1], [1, 0, 1]]
    haplotypes = solve_wMEC(allele_matrix, allow_homozygousity=True)
    assert set(haplotypes) == {(0, 1, 1), (1, 0, 1)}
    haplotypes = solve_wMEC(allele_matrix, allow_homozygousity=False)
    assert set(haplotypes) == {(0, 1, 1), (1, 0, 0)}

    allele_matrix = [[0, 1, 1], [0, 1, 1]]
    haplotypes = solve_wMEC(allele_matrix, allow_homozygousity=True)
    assert set(haplotypes) == {(0, 1, 1), (-1, -1, -1)}
    haplotypes = solve_wMEC(allele_matrix, allow_homozygousity=False)
    assert set(haplotypes) == {(0, 1, 1), (1, 0, 0)}


def test_partition():
    allele_matrix = [[0, 1, 0], [0, 1, 0], [0, 1, 1]]
    result = wMECSolver(allele_matrix).solve()
    assert result.partition == (1, 1, 0)

    allele_matrix = [[0, 1, 0], [0, 1, 0], [-1, -1, 1]]
    result = wMECSolver(allele_matrix).solve()
    assert result.partition == (1, 1, 0)

    allele_matrix = [[0, 1, 0], [0, 1, 0], [0, 1, 1]]
    result = wMECSolver(allele_matrix).solve()
    assert result.partition == (1, 1, 0)
    result = wMECSolver(allele_matrix).solve(allow_homozygousity=False)
    assert result.partition == (0, 0, 0)


def test_cost():
    allele_matrix = [[0, 1, 0], [1, 0, 1], [0, 1, 1]]
    result = wMECSolver(allele_matrix).solve()
    assert result.cost == 1
    weights = [[3, 3, 3], [3, 3, 3], [3, 3, 2.9]]  # float weights are silently rounded down
    result = wMECSolver(allele_matrix, weights).solve()
    assert result.cost == 2


def test_argsort_reads():
    reads = [(-1, 1, 0), (1, 0, -1), (-1, -1, 0)]
    assert wMECSolver._argsort_reads(reads) == [1, 0, 2]

def test_sort_reads():
    allele_matrix = [[-1, 1, 0], [1, 0, 1], [0, 1, 0]]
    weights = [[1, 1, 1], [1, 2, 1], [2, 2, 2]]
    solver = wMECSolver(allele_matrix, weights)
    assert solver._index_mapping == {0: 2, 1: 0, 2: 1}
    result = solver.solve()
    haplotypes = result.haplotypes
    assert set(haplotypes) == {(1, 0, 1), (0, 1, 0)}
    partition = result.partition
    assert partition[0] == partition[2] and partition[1] != partition[2]