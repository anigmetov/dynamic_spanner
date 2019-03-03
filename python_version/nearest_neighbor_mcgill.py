#!/usr/bin/env python3

import joblib as jl
import numpy as np
import random
import sys
inf_dist = 10000000000.0


def read_mcgill_matrix(fname):
    with open(fname, 'r') as f:
        s = f.readline();
        n_points = int(s)
        a = np.zeros((n_points, n_points))
        for i in range(n_points):
            for j in range(n_points):
                s = f.readline()
                a[i][j] = float(s)
    return (a, n_points)


class ExactSpanner:
    def __init__(self, dist_matrix):
        self.dist_matrix = dist_matrix

    def get_distance(self, p_a, p_b):
        return self.dist_matrix[p_a][p_b]

    def upper_bound_by_index(self, p_a, p_b):
        return self.dist_matrix[p_a][p_b]


class AnnFinder:
    def __init__(self, points, dyn_spanner):
        self.points = points
        self.n_points = len(points)
        self.dyn_spanner = dyn_spanner

    def discard_1(self, candidate_index):
        discard_outer_dist = self.distance_to_query(candidate_index) * (1.0 + 1.0 / (1.0 + self.epsilon))
        self.candidate_list = list(filter(
            lambda idx: self.pwdist(candidate_index, idx) < discard_outer_dist, self.candidate_list))

    def discard_2(self, worse_cand_index, candidate_index):
        dw = self.distance_to_query(worse_cand_index)
        dce = self.distance_to_query(candidate_index) / (1.0 + self.epsilon)
        discard_small_radius = dw - dce
        assert (discard_small_radius > 0.0)
        discard_large_radius = dw + dce
        discard_predicate = lambda idx: self.pwdist(worse_cand_index, idx) > discard_small_radius and self.pwdist(
            worse_cand_index, idx) < discard_large_radius
        self.candidate_list = list(filter(discard_predicate, self.candidate_list))

    def pwdist(self, idx1, idx2):
        return self.dyn_spanner.get_distance(self.points[idx1], self.points[idx2])

    def distance_to_query(self, idx):
        dist = self.dyn_spanner.get_distance(self.points[idx], self.query_point)
        self.distance_computed[idx] = True
        return dist

    def find_ann(self, query_point, epsilon):
        self.distance_computed = np.zeros(self.n_points, dtype=np.bool)
        self.query_point = query_point
        self.epsilon = epsilon
        self.candidate_list = [i for i in range(self.n_points)]
        random.shuffle(self.candidate_list)
        cand_idx = self.candidate_list.pop()
        self.discard_1(cand_idx)
        while len(self.candidate_list) > 0:
            alternative_cand_idx = self.candidate_list.pop()
            da = self.distance_to_query(alternative_cand_idx)
            dc = self.distance_to_query(cand_idx)
            if dc < da:
                self.discard_2(alternative_cand_idx, cand_idx)
            else:
                cand_idx, alternative_cand_idx = alternative_cand_idx, cand_idx
                self.discard_2(alternative_cand_idx, cand_idx)
                self.discard_1(cand_idx)

        assert(self.check_nearest_neighbor(self.points[cand_idx], epsilon))

        return self.points[cand_idx]

    def number_of_computed_distances(self):
        return np.sum(self.distance_computed)

    def fraction_of_distances_computed(self):
        return self.number_of_computed_distances() / float(self.n_points)

    def check_nearest_neighbor(self, candidate, epsilon):
        d_cand = self.dyn_spanner.get_distance(candidate, self.query_point)
        d_correct = min([self.dyn_spanner.get_distance(self.query_point, p) for p in self.points])
        return (d_cand <= (1.0 + epsilon) * d_correct)


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
            self.upper_bounds[i] = min(self.upper_bounds[i], self.dyn_spanner.upper_bound_by_index(idx, i) + dist)
        # update lower bounds
        for i in range(self.n_points):
            if self.distance_computed[i]:
                continue
            old_lb = self.lower_bounds[i]
            self.lower_bounds[i] = max(self.lower_bounds[i],
                                       dist - self.dyn_spanner.upper_bound_by_index(idx, i))
            true_dist = self.dyn_spanner.get_distance(self.query_point, self.points[i])
            # if true_dist < self.lower_bounds[i]:
            #     print(f"i = {i}, idx = {idx}, n_points = {self.n_points}, query = {self.query_point}, point[i] = {self.points[i]},  old_lb = {old_lb},  new_lb = {self.upper_bounds[i]}, distance = {true_dist}")

            # assert (self.dyn_spanner.get_distance(self.query_point, self.points[i]) >= 0.999999 * self.lower_bounds[i])

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
        assert (d_cand <= 1.0000001 * d_correct)
        # assert (d_cand == d_correct)


def run_experiment_ann(n_points, fname, epsilon, n_attempts):
    np.random.seed(1)
    random.seed(1)
    # print("reading file ", fname, "\n")
    dist_matrix, n_total_points = read_mcgill_matrix(fname)
    # print(f"n_total_points = {n_total_points}\n")
    points = [i for i in range(n_total_points)]
    np.random.shuffle(points)
    points, queries = points[:n_points], points[n_points:]
    # print(f"points = {points}; queries = {queries[0]}\n")
    spanner = ExactSpanner(dist_matrix)
    ann_finder = AnnFinder(points, spanner)
    computed_distances = np.zeros(n_attempts*n_attempts, dtype=np.int32)
    computed_distance_fractions = np.zeros(n_attempts * n_attempts, dtype=np.float)
    for j in range(n_attempts):
        query = queries[j]
        for t in range(n_attempts):
            nearest_neighbor = ann_finder.find_ann(query, epsilon)
            if not ann_finder.check_nearest_neighbor(nearest_neighbor, epsilon):
                print(f"ACHTUNG! ERROR! j = {j}, {fname};{n_points};{np.min(computed_distance_fractions)};{np.max(computed_distance_fractions)}")
            computed_distance_fractions[n_attempts * j + t] = ann_finder.fraction_of_distances_computed()
            computed_distances[n_attempts * j + t] = ann_finder.number_of_computed_distances()

    s1 = f"{fname};{epsilon};{n_points};"
    s2 = f"{np.min(computed_distance_fractions)};{np.max(computed_distance_fractions)};{np.average(computed_distance_fractions)};"
    s3 = f"{np.min(computed_distances)};{np.max(computed_distances)};{np.average(computed_distances)}"
    print(f"{s1}{s2}{s3}")
    sys.stdout.flush()
    return (n_points, fname, epsilon,
            np.min(computed_distance_fractions), np.max(computed_distance_fractions),
            np.average(computed_distance_fractions), np.min(computed_distances), np.max(computed_distances), np.average(computed_distances))



if __name__ == "__main__":
    np.random.seed(1)
    n_pointses = [100]
    # dims = [0, 1, 2]
    dims = [0]
    # qs = [1, 2, 3]
    qs = [1]
    fnames = ["dist_matrix_q_{}_dim_{}.txt".format(q, dim) for q in qs for dim in dims]
    # epsilons = [0.0001, 0.0005, 0.001, 0.005, 0.01, 0.1, 0.5]
    epsilons = [0.1]
    n_attempts = 5

    results_ann = jl.Parallel(n_jobs=-1)(
        jl.delayed(run_experiment_ann)(n_points, fname, eps, n_attempts)
        for n_points in n_pointses
        for fname in fnames
        for eps in epsilons)
    print("################################################################################")
    for r in results_ann:
        print(r)

    # results = jl.Parallel(n_jobs=-1)(
    #     jl.delayed(run_experiment)(n_points, dim, ps_gen_method, ps_arg, query_gen, n_attempts)
    #     for n_points in n_pointses
    #     for dim in dims
    #     for (ps_gen_method, ps_arg) in zip(ps_gen_methods, ps_gen_args)
    #     for query_gen in gen_queries)
    # print("################################################################################")
    # for r in results:
    #     print(r)

