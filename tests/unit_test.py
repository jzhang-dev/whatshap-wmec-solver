#!/usr/bin/env python
# coding: utf-8

from whatshap_wmec_solver import wMECSolver, solve_wMEC


def test_simple_matrix():
    # Two reads with no conflicts
    allele_matrix = [[0, 0, 1, 1], [1, 1, 0, 0]]
    assert solve_wMEC(allele_matrix) == ([0, 0, 1, 1], [1, 1, 0, 0])

    # The second variant is homozygous
    allele_matrix = [[0, 1, 0], [0, 1, 0], [0, 1, 1]]
    assert solve_wMEC(allele_matrix) == ([0, 1, 0], [0, 1, 1])

    # The last read contains an error at the second variant
    allele_matrix = [[0, 1, 0], [1, 0, 1], [1, 0, 1], [1, 1, 1]]
    assert solve_wMEC(allele_matrix) == ([0, 1, 0], [1, 0, 1])


def test_weights():
    # The last read conflicts with other reads at the second variant but has higher weights
    allele_matrix = [[0, 1, 0], [1, 0, 1], [1, 0, 1], [1, 1, 1]]
    weights = [[1, 1, 1], [1, 1, 1], [1, 1, 1], [10, 10, 10]]
    assert solve_wMEC(allele_matrix, weights=weights) == ([0, 1, 0], [1, 1, 1])


def test_ambiguous_matrix():
    # Read #3 conflicts with Read #2 at Variant #2. Without additional evidence, Variant #2 is ambiguous. 
    allele_matrix = [[0, 1, 0], [1, 0, 1], [1, 1, 1]]
    assert solve_wMEC(allele_matrix) == ([0, 1, 0], [1, -1, 1])

def test_reads_with_holes():
    allele_matrix = [[0, 1, 0], [1, 0, 1], [1, -1, 1]]
    assert solve_wMEC(allele_matrix) == ([0, 1, 0], [1, 0, 1])


def test_multiple_phase_blocks():
    allele_matrix = [[0, 1, -1, -1], [1, 0, -1, -1], [-1, -1, 1, 0], [-1, -1, 0, 1]]
    assert solve_wMEC(allele_matrix) == ([0, 1, 1, 0], [1, 0, 0, 1])


def test_allow_homozygous():
    allele_matrix = [[0, 1, 1], [0, 1, 1], [1, 0, 1]]
    assert solve_wMEC(allele_matrix, allow_homozygousity=True) == ([0, 1, 1], [1, 0, 1])
    assert solve_wMEC(allele_matrix, allow_homozygousity=False) == ([0, 1, 1], [1, 0, 0])

    allele_matrix = [[0, 1, 1], [0, 1, 1]]
    assert solve_wMEC(allele_matrix, allow_homozygousity=True) == ([-1, -1, -1], [0, 1, 1])
    assert solve_wMEC(allele_matrix, allow_homozygousity=False) == ([0, 1, 1], [1, 0, 0])
