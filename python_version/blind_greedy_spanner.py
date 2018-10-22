#!/usr/bin/env python3

import networkx as nx
import math
import numpy as np
import random
import sys

from utils_euclidean import *

inf_dist = 10000000000.0


class BlindGreedySpanner:

    def __init__(self, points, inf_dist, dist_function=None, dist_matrix=None, compute_distance_matrix=False):
        assert (dist_function != None or dist_matrix != None)
        self.points = points
        self.n_points = len(points)

        # common shape for all matrices
        matr_shape = (n_points, n_points)

        if dist_matrix != None:
            self.dist_matrix = dist_matrix
        else:
            self.dist_matrix = np.zeros(matr_shape)
            if compute_distance_matrix:
                for i in range(self.n_points):
                    for j in range(i + 1, self.n_points):
                        self.dist_matrix[i][j] = dist_function(points[i], points[j])

        self.lower_bounds = np.zeros(matr_shape)
        self.upper_bounds = np.full(matr_shape, inf_dist)

        self.dist_requested = np.zeros(matr_shape, dtype=np.int32)
        self.dist_computed = np.zeros(matr_shape, dtype=np.int32)
        self.dist_function = dist_function

        self.spanner = nx.Graph()
        for i in range(self.n_points):
            self.spanner.add_node(i)


    def set_upper_bound(self, i, j, value):
        assert(value >= 0.0 and value <= self.upper_bounds[i][j])
        self.upper_bounds[i][j] = self.upper_bounds[j][i] = value

    def set_lower_bound(self, i, j, value):
        assert(value >= self.lower_bounds[i][j])
        self.lower_bounds[i][j] = self.lower_bounds[j][i] = value

    def upper_bound(self, i, j):
        return self.upper_bounds[i][j]

    def distance(self, i, j):
        if i == j:
            return 0.0
        self.dist_requested[i][j] = self.dist_requested[j][i] = 1
        self.dist_computed[i][j] = self.dist_computed[j][i] = 1
        result = self.dist_matrix[i][j]
        if result == 0.0:
            result = self.dist_function(self.points[i], self.points[j])
            self.dist_matrix[i][j] = self.dist_matrix[j][i] = result
        self.set_lower_bound(i, j, result)
        self.set_upper_bound(i, j, result)
        return result


    def worst_pair_blind(self):
        i = 0
        j = 1
        worst_ratio = inf_dist
        return (worst_ratio, i, j)

    def add_edge_to_spanner(self, i, j):
        dist_ij = self.distance(i, j)
        assert(dist_ij > 0.0 and i != j)
        self.spanner.add_weighted_edges_from([i, j, dist_ij])
        self.update_bounds_using_distance(i, j)


    def update_bounds_using_distance(self, i, j):
        self.update_upper_bounds(i, j)
        self.update_lower_bounds_1(i, j)
        self.update_lower_bounds_2(i, j)
        # self.update_lower_bounds_1(i, j)
        # self.update_lower_bounds_2(i, j)


    def update_upper_bounds(self, i, j):
        dist_ij = self.dist_matrix[i][j]
        assert(dist_ij > 0.0)

        for x in range(self.n_points):
            for y in range(x + 1, self.n_points):
                if self.dist_computed[x][y]:
                    continue
                new_ub = min(self.upper_bounds[x][y],
                             self.upper_bounds[x][i] + dist_ij + self.upper_bounds[y][j],
                             self.upper_bounds[x][j] + dist_ij + self.upper_bounds[y][i])
                self.set_upper_bound(i, j, new_ub)


    def update_lower_bounds_1(self, i, j):
        dist_ij = self.dist_matrix[i][j]
        assert(dist_ij > 0.0)

        for x in range(self.n_points):
            for y in range(x + 1, self.n_points):
                if self.dist_computed[x][y]:
                    continue
                new_lb = max(self.lower_bounds[i][j],
                           dist_ij - self.upper_bounds[x][i] - self.upper_bounds[y][j],
                           dist_ij - self.upper_bounds[x][j] - self.upper_bounds[y][i])
                self.set_lower_bound(i, j, new_lb)


    def update_lower_bounds_2(self, i, j):
        dist_ij = self.dist_matrix[i][j]
        assert(dist_ij > 0.0)
        for x in range(self.n_points):
            for y in range(x + 1, self.n_points):
                if self.dist_computed[x][y]:
                    continue
                new_lb = max(self.lower_bounds[i][j],
                             self.lower_bounds[j][y] - dist_ij - self.upper_bounds[x][i],
                             self.lower_bounds[x][i] - dist_ij - self.upper_bounds[y][j],
                             self.lower_bounds[i][y] - dist_ij - self.upper_bounds[x][j],
                             self.lower_bounds[x][j] - dist_ij - self.upper_bounds[y][i])
                self.set_lower_bound(i, j, new_lb)

    def build_blind(self, epsilon):
        while True:
            worst_ratio, i, j = self.worst_pair_blind()
            if worst_ratio <= (1 + epsilon):
                break
            self.add_edge_to_spanner(i, j)



def distance_matrix(points):
    n_points = len(points)
    result = np.zeros((n_points, n_points))
    for i in range(n_points):
        for j in range(n_points):
            if i <= j:
                continue
            result[i][j] = result[j][i] = euclidean_distance(points[i, :], points[j, :])
    return result


# def distances_list(points):
#     n_points = len(points)
#     dm = distance_matrix(points)
#     return [(dm[i][j], i, j) for i in range(n_points) for j in range(i + 1, n_points)]
#
#
# def sorted_distances(points):
#     return sorted(distances_list(points))
#
#
# def upper_bound(g, i, j):
#     try:
#         return nx.shortest_path_length(g, i, j, weight='weight')
#     except nx.NetworkXNoPath:
#         return inf_dist
#
#
# def get_distance(a, b):
#     return np.linalg.norm(a - b)


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


# def get_buckets(max_dist, n_buckets):
#     scaling_factor = max_dist / (2 ** (n_buckets))
#     result = [(0.0, scaling_factor * 2.0)]
#     result.extend([(scaling_factor * (2 ** k), 4.0 * scaling_factor * (2 ** k)) for k in range(n_buckets - 1)])
#     return result
#
#
# def get_bucket_index(buckets, d):
#     dist = d[0]
#     cand_list = []
#     for i in range(len(buckets)):
#         lb, ub = buckets[i]
#         if dist >= lb and dist <= ub:
#             cand_list.append(i)
#     if len(cand_list) == 0:
#         b0 = buckets[0]
#         b1 = buckets[-1]
#         print(f"ACHTUNG! dist = {dist}, {b0} - {b1}  bucket not found")
#     assert (len(cand_list) > 0)
#     if len(cand_list) == 1:
#         return cand_list[0]
#     else:
#         idx = random.randint(0, len(cand_list) - 1)
#         return cand_list[idx]
#
#
# def quasi_sorted_distances(points, n_buckets, dist_list=None):
#     if dist_list is None:
#         n_points = len(points)
#         dm = distance_matrix(points)
#         dist_list = [(dm[i][j], i, j) for i in range(n_points) for j in range(i + 1, n_points)]
#     max_dist = max(dist_list)[0]
#     buckets = get_buckets(max_dist, n_buckets)
#     result = [[] for i in range(n_buckets)]
#     for d in dist_list:
#         result[get_bucket_index(buckets, d)].append(d)
#     map(random.shuffle, result)
#     return [x for y in result for x in y]
#
#
# def blind_greedy_spanner(points, eps):
#     return
#
#
# def greedy_spanner(points, eps, sd=None):
#     t = 1.0 + eps
#     if sd is None:
#         sd = sorted_distances(points)
#     spanner = nx.Graph()
#     for i in range(len(points)):
#         spanner.add_node(i)
#     for (d, i, j) in sd:
#         ub = upper_bound(spanner, i, j)
#         if ub > t * d:
#             spanner.add_edge(i, j, weight=d)
#     return spanner


if __name__ == "__main__":

    # We like goblinses, batses and fishes. But we hasn't tried Hobbitses before.
    n_pointses = [1000]
    # n_bucketses = [5, 10, 15]
    epsilons = [0.1]
    dims = [40, 80]

    for dim in dims:
        for n_points in n_pointses:
            points = get_points(n_points, dim)
            dist_list = distances_list(points)
            min_dist = min(dist_list)[0]
            max_dist = max(dist_list)[0]
            spread = max_dist / min_dist
            n_buckets = int(math.log(spread, 2.0))
            print(
                f"dim = {dim}, #points = {n_points}, spread = {max_dist} / {min_dist} =  {spread}, n_buckets = {n_buckets}")
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

                sd = quasi_sorted_distances(points, n_buckets, dist_list)
                gs = greedy_spanner(points, eps, sd)
                n_gs_quasi_edges = gs.number_of_edges()
                sparseness_quasi = n_gs_quasi_edges / total_n_pairs
                print(
                    f"dim = {dim}, eps = {eps}, #points = {n_points}, #buckets =  {n_buckets} with quasi-sorted distances sparseness = {n_gs_quasi_edges} / {total_n_pairs}  = {sparseness_quasi}")

    sys.exit(0)
