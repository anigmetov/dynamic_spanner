#!/usr/bin/env python3

import networkx as nx
import numpy as np
import random
import sys

inf_dist = 10000000000


class QuasiSortedDistances:
    def __init__(self, qsorted_distances, n_buckets):
        self.qsorted_distances = qsorted_distances
        self.n_buckets = n_buckets

def get_points(n_points, dim):
    return np.random.randn(n_points, dim)


def get_exponential_points(n_points, dim):
    max_power = max(n_points / 2.0, 10.0)
    print(f"max_power = {max_power}")
    shape = (n_points, dim)
    res = np.random.uniform(1.0, max_power, shape)
    return np.power(2.0, res)


def distance_matrix(points):
    n_points = len(points)
    result = np.zeros((n_points, n_points))
    for i in range(n_points):
        for j in range(n_points):
            if i <= j:
                continue
            result[i][j] = result[j][i] = get_distance(points[i, :], points[j, :])
    return result


def distances_list(points):
    n_points = len(points)
    dm = distance_matrix(points)
    return [(dm[i][j], i, j) for i in range(n_points) for j in range(i + 1, n_points)]


def sorted_distances(points):
    return sorted(distances_list(points))


def upper_bound(g, i, j):
    try:
        return nx.shortest_path_length(g, i, j, weight='weight')
    except nx.NetworkXNoPath:
        return inf_dist


def get_distance(a, b):
    return np.linalg.norm(a - b)


def get_buckets(min_dist, max_dist):
    n_buckets = int(np.ceil(np.log2(max_dist/min_dist)))
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


if __name__ == "__main__":

    # min_dist = 0.01
    # max_dist = 3.0
    # print(get_buckets(min_dist, max_dist))

    n_pointses = [50, 100]
    epsilons = [0.01, 0.1, 0.5]
    dims = [3, 5, 7]

    for dim in dims:
        for n_points in n_pointses:
            # points = get_points(n_points, dim)
            points = get_exponential_points(n_points, dim)
            dist_list = distances_list(points)
            max_dist = max(dist_list)[0]
            for eps in epsilons:

                gs = greedy_spanner(points, eps, sorted(dist_list))
                total_n_pairs = n_points * (n_points - 1) / 2.0
                n_gs_edges = gs.number_of_edges()
                sparseness = n_gs_edges / total_n_pairs
                print(
                    f"dim = {dim}, eps = {eps}, #points = {n_points}, greedy spanner sparseness = {n_gs_edges} / {total_n_pairs}  = {sparseness}")

                random_dist_list = dist_list[:]
                random.shuffle(random_dist_list)
                gs_random = greedy_spanner(points, eps, random_dist_list)

                n_gs_random_edges = gs_random.number_of_edges()
                sparseness_random = n_gs_random_edges / total_n_pairs
                print(
                    f"dim = {dim}, eps = {eps}, #points = {n_points}, random order spanner sparseness = {n_gs_random_edges} / {total_n_pairs}  = {sparseness_random}")

                sd = quasi_sorted_distances(points, dist_list)
                gs = greedy_spanner(points, eps, sd.qsorted_distances)
                n_gs_quasi_edges = gs.number_of_edges()
                sparseness = n_gs_quasi_edges / total_n_pairs
                print(
                    f"dim = {dim}, eps = {eps}, #points = {n_points}, #buckets =  {sd.n_buckets} with quasi-sorted distances sparseness = {n_gs_quasi_edges} / {total_n_pairs}  = {sparseness}")

    sys.exit(0)