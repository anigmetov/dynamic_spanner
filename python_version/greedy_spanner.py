#!/usr/bin/env python3

import networkx as nx
import numpy as np
import random
import sys

inf_dist = 10000000000

dim = 2


def get_points(n_points):
    return np.random.randn(n_points, dim)


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


# def is_distance_larger(a, b, t):
#     """Return true, if d(a, b) > t;
#        return false, if d(a,b) < t/2; if d(a,b) is between t/2 and t,
#        return either true or false"""
#     d = get_distance(a, b)
#     if d > t:
#         return True
#     elif d < t / 2:
#         return False
#     else:
#         return np.random.uniform() < 0.5


def get_buckets(max_dist, n_buckets):
    scaling_factor = max_dist / (2 ** (n_buckets))
    result = [(0.0, scaling_factor * 2.0)]
    result.extend([(scaling_factor * (2 ** k), 4.0 * scaling_factor * (2 ** k)) for k in range(n_buckets - 1)])
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


def quasi_sorted_distances(points, n_buckets, dist_list=None):
    if dist_list is None:
        n_points = len(points)
        dm = distance_matrix(points)
        dist_list = [(dm[i][j], i, j) for i in range(n_points) for j in range(i + 1, n_points)]
    max_dist = max(dist_list)[0]
    buckets = get_buckets(max_dist, n_buckets)
    result = [[] for i in range(n_buckets)]
    for d in dist_list:
        result[get_bucket_index(buckets, d)].append(d)
    map(random.shuffle, result)
    return [x for y in result for x in y]


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

    # We like goblinses, batses and fishes. But we hasn't tried Hobbitses before.
    n_pointses = [50, 100, 200, 500, 1000, 2000]
    n_bucketses = [10, 20, 30]
    epsilons = [0.01, 0.1, 0.2, 0.5]

    for n_points in n_pointses:
        points = get_points(n_points)
        dist_list = distances_list(points)
        for eps in epsilons:

            gs = greedy_spanner(points, eps, sorted(dist_list))
            total_n_pairs = n_points * (n_points - 1) / 2.0
            n_gs_edges = gs.number_of_edges()
            sparseness = n_gs_edges / total_n_pairs
            print(
                f"eps = {eps}, #points = {n_points}, greedy spanner sparseness = {n_gs_edges} / {total_n_pairs}  = {sparseness}")

            for n_buckets in n_bucketses:
                sd = quasi_sorted_distances(points, n_buckets, dist_list)
                gs = greedy_spanner(points, eps, sd)
                n_gs_quasi_edges = gs.number_of_edges()
                sparseness = n_gs_quasi_edges / total_n_pairs
                print(
                    f"eps = {eps}, #points = {n_points}, #buckets =  {n_buckets} with quasi-sorted distances sparseness = {n_gs_quasi_edges} / {total_n_pairs}  = {sparseness}")

    sys.exit(0)
