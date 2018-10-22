#!/usr/bin/env python3

import joblib as jl
import numpy as np
import random
from utils_euclidean import *
from blind_greedy_spanner import *

inf_dist = 10000000000.0


class ExactSpanner:
    def __init__(self, points):
        self.points = points
        self.n_points = len(points)
        self.dist_matrix = distance_matrix(points)

    def get_distance(self, p_a, p_b):
        return euclidean_distance(p_a, p_b)

    def upper_bound_by_index(self, i, j):
        return euclidean_distance(self.points[i], self.points[j])


class NearestNeighborFinder:
    def __init__(self, points, dyn_spanner):
        self.points = points
        self.n_points = len(points)
        self.dyn_spanner = dyn_spanner
        self.reset_bounds()

    def reset_bounds(self):
        self.upper_bounds = np.full(self.n_points, inf_dist)
        self.lower_bounds = np.zeros(self.n_points)
        self.distances = np.full(self.n_points, -1.0)
        self.distance_computed = np.zeros(self.n_points)

    def smallest_upper_bound(self):
        return np.argmin(self.upper_bounds)

    def get_better_candidates(self, candidate_index):
        return np.nonzero(self.lower_bounds < self.upper_bounds[candidate_index])[0]

    def compute_distance(self, idx):
        dist = self.dyn_spanner.get_distance(self.points[idx], self.query_point)
        self.distances[idx] = dist
        self.distance_computed[idx] = 1
        self.upper_bounds[idx] = dist
        self.lower_bounds[idx] = dist
        self.update_bounds(idx)
        return dist

    def update_bounds(self, idx):
        dist = self.distances[idx]
        # update upper bounds
        for i in range(self.n_points):
            if self.distance_computed[i]:
                continue
            ub = self.dyn_spanner.upper_bound_by_index(idx, i)
            assert (ub == self.dyn_spanner.get_distance(self.points[i], self.points[idx]))
            assert (dist == self.dyn_spanner.get_distance(self.query_point, self.points[idx]))
            self.upper_bounds[i] = min(self.upper_bounds[i], self.dyn_spanner.upper_bound_by_index(idx, i) + dist)
            # true_dist = self.dyn_spanner.get_distance(self.query_point, self.points[i])
            # if true_dist > self.upper_bounds[i]:
            #     print(f"query = {self.query_point}, point[i] = {self.points[i]},  new ub = {self.upper_bounds[i]}, distance = {true_dist}")
            assert (self.dyn_spanner.get_distance(self.query_point, self.points[i]) <= self.upper_bounds[i])

            # update lower bounds
        for i in range(self.n_points):
            if self.distance_computed[i]:
                continue
            self.lower_bounds[i] = max(self.lower_bounds[i],
                                       dist - self.dyn_spanner.upper_bound_by_index(idx, i))
            assert (self.dyn_spanner.get_distance(self.query_point, self.points[i]) >= self.lower_bounds[i])

    def find_nn(self, query_point):
        self.reset_bounds()
        self.query_point = query_point
        while True:
            current_candidate_index = self.smallest_upper_bound()
            better_candidates = self.get_better_candidates(current_candidate_index)
            if better_candidates.size:
                self.compute_distance(random.choice(better_candidates))
            else:
                return self.points[current_candidate_index]

    def number_of_computed_distances(self):
        return np.sum(self.distance_computed)

    def fraction_of_distances_computed(self):
        return self.number_of_computed_distances() / float(self.n_points)

    def check_nearest_neighbor(self, candidate):
        d_cand = self.dyn_spanner.get_distance(candidate, self.query_point)
        d_correct = min([self.dyn_spanner.get_distance(self.query_point, p) for p in self.points])
        # assert (d_cand <= 1.0000001 * d_correct)
        assert (d_cand == d_correct)


def run_experiment(n_points, dim, pointset_generation_method, args, query_generation_method, n_attempts):
    points = pointset_generation_method(n_points, dim, *args)
    dist_matrix = distance_matrix(points)
    # spanner = BlindGreedySpanner(points, inf_dist, dist_matrix=dist_matrix, dist_function=euclidean_distance)
    spanner = ExactSpanner(points)
    nn_finder = NearestNeighborFinder(points, spanner)
    computed_distance_fractions = np.zeros(n_attempts, dtype=np.float)
    for j in range(n_attempts):
        query = query_generation_method(dim)
        nearest_neighbor = nn_finder.find_nn(query)
        nn_finder.check_nearest_neighbor(nearest_neighbor)
        computed_distance_fractions[j] = nn_finder.fraction_of_distances_computed()

    print(f"{dim};{n_points};{pointset_generation_method.__name__};{np.min(computed_distance_fractions)};{np.max(computed_distance_fractions)}")
    return (n_points, dim, pointset_generation_method.__name__, args, query_generation_method.__name__,
            np.min(computed_distance_fractions), np.max(computed_distance_fractions), np.average(computed_distance_fractions))


def query_1(dim):
    return np.random.uniform(-1000.0, 1000.0, (1, dim))


def query_2(dim):
    return np.random.normal(scale=100.0, size=(1, dim))


if __name__ == "__main__":
    np.random.seed(1)
    # n_pointses = [50]
    n_pointses = [50, 100, 500, 1000, 5000, 10000]
    # dims = [2]
    dims = [2,3,10,20]
    ps_gen_methods = [get_points, get_uniform_points, get_exponential_points]
    ps_gen_args = [[], [10.0], []]
    gen_queries = [query_1, query_2]
    n_attempts = 10

    results = jl.Parallel(n_jobs=-1)(jl.delayed(run_experiment)(n_points, dim, ps_gen_method, ps_arg, query_gen, n_attempts)
                                    for n_points in n_pointses
                                    for dim in dims
                                    for (ps_gen_method, ps_arg) in zip(ps_gen_methods, ps_gen_args)
                                    for query_gen in gen_queries)
    print("################################################################################")
    for r in results:
        print(r)
