#!/usr/bin/env python3

import networkx as nx
import numpy as np
import random
import sys
import matplotlib.pyplot as plt
import joblib as jl

from utils_euclidean import *

inf_dist = 10000000000


class QuasiSortedDistances:
    def __init__(self, qsorted_distances, n_buckets):
        self.qsorted_distances = qsorted_distances
        self.n_buckets = n_buckets


def upper_bound(g, i, j):
    try:
        return nx.shortest_path_length(g, i, j, weight='weight')
    except nx.NetworkXNoPath:
        return inf_dist


def get_buckets(min_dist, max_dist):
    n_buckets = int(np.ceil(np.log2(max_dist / min_dist)))
    scaling_factor = min_dist
    result = [(scaling_factor * (2 ** k), 4.0 * scaling_factor * (2 ** k)) for k in range(n_buckets)]
    return result


def get_bucket_index(buckets, d):
    dist = d[0]
    cand_list = []
    for i in range(len(buckets)):
        lb, ub = buckets[i]
        if dist >= lb and dist <= ub:
            cand_list.append(i)
    if len(cand_list) == 0:
        b0 = buckets[0]
        b1 = buckets[-1]
        print(f"ACHTUNG! dist = {dist}, {b0} - {b1}  bucket not found")
    assert (len(cand_list) > 0)
    if len(cand_list) == 1:
        return cand_list[0]
    else:
        idx = random.randint(0, len(cand_list) - 1)
        return cand_list[idx]


def quasi_sorted_distances(points, dist_list=None):
    if dist_list is None:
        n_points = len(points)
        dm = distance_matrix(points)
        dist_list = [(dm[i][j], i, j) for i in range(n_points) for j in range(i + 1, n_points)]
    max_dist = max(dist_list)[0]
    min_dist = min(dist_list)[0]
    buckets = get_buckets(min_dist, max_dist)
    n_buckets = len(buckets)
    result = [[] for i in range(n_buckets)]
    for d in dist_list:
        result[get_bucket_index(buckets, d)].append(d)
    map(random.shuffle, result)
    return QuasiSortedDistances([x for y in result for x in y], n_buckets)


def greedy_spanner(points, eps, sd=None):
    t = 1.0 + eps
    if sd is None:
        sd = sorted_distances(points)
    spanner = nx.Graph()
    for i in range(len(points)):
        spanner.add_node(i)
    for (d, i, j) in sd:
        ub = upper_bound(spanner, i, j)
        if ub > t * d:
            spanner.add_edge(i, j, weight=d)
    return spanner


class ExperimentResult:
    def __init__(self, dim, n_points, eps, points_gen_method, spanner_method, spanner_edges):
        self.dim = dim
        self.n_points = n_points
        self.eps = eps
        self.total_n_pairs = n_points * (n_points - 1) / 2.0
        self.spanner_method = spanner_method
        self.spanner_edges = spanner_edges
        self.sparseness = spanner_edges / self.total_n_pairs
        self.points_method = points_gen_method

    def __str__(self):
        return f"dim = {self.dim}, eps = {self.eps}, #points = {self.n_points}, points method = {self.points_method}, spanner = {self.spanner_method}, sparseness = {self.spanner_edges} / {self.total_n_pairs}  = {self.sparseness}"


def run_experiment(dim, n_points, points_generator, eps):
    points = points_generator(n_points, dim)
    dist_list = distances_list(points)

    result = []

    gs = greedy_spanner(points, eps, sorted(dist_list))
    n_gs_edges = gs.number_of_edges()

    result.append(ExperimentResult(dim, n_points, eps, points_generator.__name__, "greedy", n_gs_edges))

    print(result[-1])

    random_dist_list = dist_list[:]
    random.shuffle(random_dist_list)
    gs_random = greedy_spanner(points, eps, random_dist_list)

    n_gs_random_edges = gs_random.number_of_edges()
    result.append(ExperimentResult(dim, n_points, eps, points_generator.__name__, "random", n_gs_random_edges))
    print(result[-1])

    sd = quasi_sorted_distances(points, dist_list)
    gs_quasi = greedy_spanner(points, eps, sd.qsorted_distances)
    n_gs_quasi_edges = gs_quasi.number_of_edges()
    result.append(ExperimentResult(dim, n_points, eps, points_generator.__name__, "quasi", n_gs_quasi_edges))
    print(result[-1])

    sys.stdout.flush()

    return result


if __name__ == "__main__":

    # min_dist = 0.01
    # max_dist = 3.0
    # print(get_buckets(min_dist, max_dist))

    n_pointses = [100]
    epsilons = [0.1]
    dims = [2]
    points_generators = [get_points, get_exponential_points, get_uniform_points]

    all_results = jl.Parallel(n_jobs=5)(jl.delayed(run_experiment)(dim, n_points, points_generator, eps)
                                         for dim in dims
                                         for n_points in n_pointses
                                         for points_generator in points_generators
                                         for eps in epsilons)

    print("##################################################")

    for r in all_results:
        for er in r:
            print(er)


    sys.exit(0)
