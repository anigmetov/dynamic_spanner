#!/usr/bin/env python3

import networkx as nx
import math
import numpy as np
import random
import sys
import joblib as jl
import pandas as pd

from spanner_experiment_result import ExperimentResult
from utils_euclidean import *
from weakly_sorted import quasi_sorted_distances

inf_dist = 10000000000.0

class BlindSpanner:

    def __init__(self, points, inf_dist, dist_function=None, dist_matrix=None, compute_distance_matrix=False):
        assert (dist_function != None or dist_matrix != None)
        self.points = points
        self.n_points = len(points)

        # common shape for all matrices
        self.clear_spanner()
        matr_shape = (self.n_points, self.n_points)

        if dist_matrix != None:
            self.dist_matrix = dist_matrix
        else:
            self.dist_matrix = np.zeros(matr_shape)
            if compute_distance_matrix:
                for i in range(self.n_points):
                    for j in range(i + 1, self.n_points):
                        self.dist_matrix[i][j] = dist_function(points[i], points[j])

        self.dist_function = dist_function

    def clear_spanner(self):
        # common shape for all matrices
        matr_shape = (self.n_points, self.n_points)

        self.lower_bounds = np.zeros(matr_shape)
        self.upper_bounds = np.full(matr_shape, inf_dist)

        self.ratios = self.upper_bounds / self.lower_bounds
        np.fill_diagonal(self.ratios, -np.inf)

        self.worst_ratio = np.inf
        self.worst_ratio

        self.dist_requested = np.zeros(matr_shape, dtype=np.bool)
        self.dist_computed = np.eye(self.n_points, dtype=np.bool)

        self.spanner = nx.Graph()
        for i in range(self.n_points):
            self.spanner.add_node(i)

    def set_upper_bound(self, i, j, value):
        assert (value >= 0.0 and value <= self.upper_bounds[i][j])
        # if value < self.get_distance_no_cache(i, j):
        #     print(f"i = {i}, j = {j}, value = {value}, true distance = {self.get_distance_no_cache(i, j)}")
        assert (value >= self.get_distance_no_cache(i, j))
        self.upper_bounds[i][j] = self.upper_bounds[j][i] = value

    def set_lower_bound(self, i, j, value):
        assert (value >= self.lower_bounds[i][j])
        assert (value <= self.get_distance_no_cache(i, j))
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
        # print(f"in distance, result = {result}, i = {i}, j = {j}")
        self.set_lower_bound(i, j, result)
        self.set_upper_bound(i, j, result)
        return result

    def get_distance_no_cache(self, i, j):
        result = self.dist_matrix[i][j]
        if result == 0.0:
            result = self.dist_function(self.points[i], self.points[j])
        return result


    def worst_pair_blind(self):
        np.seterr(divide='ignore')
        ratio_matr = self.upper_bounds / self.lower_bounds
        # ratio_matr[ratio_matr == np.inf] = -np.inf
        np.fill_diagonal(ratio_matr, -np.inf)
        worst_ratio = ratio_matr.max()
        row, col = np.where(ratio_matr == worst_ratio)
        i = np.random.randint(len(row))
        return (worst_ratio, row[i], col[i])


    def is_blind_spanner_done(self, epsilon):
        np.seterr(divide='ignore')
        ratio_matr = self.upper_bounds / self.lower_bounds
        np.fill_diagonal(ratio_matr, -np.inf)
        return (ratio_matr <= epsilon).all()

    def add_edge_to_spanner(self, i, j):
        dist_ij = self.distance(i, j)
        assert (dist_ij > 0.0 and i != j)
        self.spanner.add_weighted_edges_from([(i, j, dist_ij)])
        self.update_bounds_using_distance(i, j)

    def number_of_edges(self):
        return self.spanner.number_of_edges()

    def update_bounds_using_distance(self, i, j):
        # old_upper_bounds = np.copy(self.upper_bounds)
        # old_lower_bounds = np.copy(self.lower_bounds)

        self.update_upper_bounds(i, j)
        self.update_lower_bounds_1(i, j)
        self.update_lower_bounds_2(i, j)

        # new_upper_bounds_1 = np.copy(self.upper_bounds)
        # new_lower_bounds_1 = np.copy(self.lower_bounds)

        # self.update_lower_bounds_1(i, j)
        # self.update_lower_bounds_2(i, j)

        # print(f"Upper bounds do not change after second run: {(self.upper_bounds == new_upper_bounds_1).all()}")
        # print(f"Lower bounds do not change after second run: {(self.lower_bounds == new_lower_bounds_1).all()}")

    def update_upper_bounds(self, i, j):
        dist_ij = self.dist_matrix[i][j]
        assert (dist_ij > 0.0)
        assert (dist_ij < inf_dist)

        for x in range(self.n_points):
            for y in range(x + 1, self.n_points):
                if self.dist_computed[x][y]:
                    continue
                new_ub = min(self.upper_bounds[x][y],
                             self.upper_bounds[x][i] + dist_ij + self.upper_bounds[y][j],
                             self.upper_bounds[x][j] + dist_ij + self.upper_bounds[y][i])
                self.set_upper_bound(x, y, new_ub)

    def update_lower_bounds_1(self, i, j):
        dist_ij = self.dist_matrix[i][j]
        assert (dist_ij > 0.0)

        for x in range(self.n_points):
            for y in range(x + 1, self.n_points):
                if self.dist_computed[x][y]:
                    continue
                new_lb = max(self.lower_bounds[x][y],
                             dist_ij - self.upper_bounds[x][i] - self.upper_bounds[y][j],
                             dist_ij - self.upper_bounds[x][j] - self.upper_bounds[y][i])
                self.set_lower_bound(x, y, new_lb)

    def update_lower_bounds_2(self, i, j):
        dist_ij = self.dist_matrix[i][j]
        assert (dist_ij > 0.0)
        for x in range(self.n_points):
            for y in range(x + 1, self.n_points):
                if self.dist_computed[x][y]:
                    continue
                new_lb = max(self.lower_bounds[x][y],
                             self.lower_bounds[j][y] - dist_ij - self.upper_bounds[x][i],
                             self.lower_bounds[x][i] - dist_ij - self.upper_bounds[y][j],
                             self.lower_bounds[i][y] - dist_ij - self.upper_bounds[x][j],
                             self.lower_bounds[x][j] - dist_ij - self.upper_bounds[y][i])
                self.set_lower_bound(x, y, new_lb)

    def build_greedy(self, epsilon):
        self.clear_spanner()
        while True:
            worst_ratio, i, j = self.worst_pair_blind()
            if worst_ratio <= (1 + epsilon):
                break
            self.add_edge_to_spanner(i, j)

    def build_random(self, epsilon):
        self.clear_spanner()
        edge_list = [(i, j) for i in range(self.n_points) for j in range(i+1, self.n_points)]
        random_edge_list = np.random.permutation(edge_list).tolist()
        while True:
            worst_ratio, _, _ = self.worst_pair_blind()
            if worst_ratio <= (1 + epsilon):
                break
            random_edge = random_edge_list.pop()
            self.add_edge_to_spanner(random_edge[0], random_edge[1])


    def build_quasi_sorted_greedy(self, epsilon):
        self.clear_spanner()
        dist_list = distances_list(self.points)
        sd = quasi_sorted_distances(self.points, dist_list)
        for dist_ij, i, j in sd.qsorted_distances:
            if self.is_blind_spanner_done(epsilon):
                break
            self.add_edge_to_spanner(i,j)


    def build_quasi_sorted_shaker(self, epsilon):
        self.clear_spanner()
        dist_list = distances_list(self.points)
        sd = quasi_sorted_distances(self.points, dist_list).qsorted_distances
        idx_min = 0
        idx_max = len(sd) - 1
        go_up = True
        for dist_ij, i, j in sd:
            if self.is_blind_spanner_done(epsilon):
                break
            if go_up:
                _, i, j = sd[idx_min]
                idx_min += 1
            else:
                _, i, j = sd[idx_max]
                idx_max -= 1
            go_up = not go_up
            self.add_edge_to_spanner(i,j)


def distance_matrix(points):
    n_points = len(points)
    result = np.zeros((n_points, n_points))
    for i in range(n_points):
        for j in range(n_points):
            if i <= j:
                continue
            result[i][j] = result[j][i] = get_lp_distance(points[i, :], points[j, :])
    return result


def run_experiment(dim, n_points, epsilon, points_generator, ps_gen_args):
    points = points_generator(n_points, dim, *ps_gen_args)
    result = []
    gs = BlindSpanner(points, inf_dist, dist_function=get_lp_distance, dist_matrix=None)
    gs.build_random(epsilon)
    n_gs_edges = gs.number_of_edges()
    result.append(ExperimentResult(dim, n_points, epsilon, points_generator.__name__, "blind_random", n_gs_edges))
    print(result[-1])
    sys.stdout.flush()
    return result

def run_experiment_quasi_greedy(dim, n_points, epsilon, points_generator, ps_gen_args):
    points = points_generator(n_points, dim, *ps_gen_args)
    result = []
    gs = BlindSpanner(points, inf_dist, dist_function=get_lp_distance, dist_matrix=None)
    gs.build_quasi_sorted_greedy(epsilon)
    n_gs_edges = gs.number_of_edges()
    result.append(ExperimentResult(dim, n_points, epsilon, points_generator.__name__, "blind_quasi_sorted_greedy", n_gs_edges))
    print(result[-1])
    sys.stdout.flush()
    gs.clear_spanner()
    gs.build_quasi_sorted_shaker(epsilon)
    n_gs_edges = gs.number_of_edges()
    result.append(ExperimentResult(dim, n_points, epsilon, points_generator.__name__, "blind_quasi_sorted_shaker", n_gs_edges))
    print(result[-1])
    sys.stdout.flush()
    return result

def run_experiment_quasi_shaker(dim, n_points, epsilon, points_generator, ps_gen_args):
    points = points_generator(n_points, dim, *ps_gen_args)
    result = []
    gs = BlindSpanner(points, inf_dist, dist_function=get_lp_distance, dist_matrix=None)
    gs.build_quasi_sorted_shaker(epsilon)
    n_gs_edges = gs.number_of_edges()
    result.append(ExperimentResult(dim, n_points, epsilon, points_generator.__name__, "blind_quasi_sorted_shaker", n_gs_edges))
    print(result[-1])
    sys.stdout.flush()
    return result

if __name__ == "__main__":
    np.random.seed(1)

    # n_pointses = [70]
    # dims = [2]
    # ps_gen_methods = [get_points, get_uniform_points]
    # ps_gen_args = [[], [10.0]]
    # epsilons = [0.2]

    n_pointses = [100, 200, 400, 600]
    dims = [2, 3]
    ps_gen_methods = [get_points, get_uniform_points]
    ps_gen_args = [[], [10.0]]
    epsilons = [0.1, 0.5]
    # run_experiment_quasi_greedy(2, 20, 0.5, get_points, [])

    results = jl.Parallel(n_jobs=-1)(
        jl.delayed(run_experiment_quasi_greedy)(dim, n_points, epsilon, ps_gen_method, ps_arg)
        for dim in dims
        for n_points in n_pointses
        for epsilon in epsilons
        for (ps_gen_method, ps_arg) in zip(ps_gen_methods, ps_gen_args))

    print("################################################################################")

    df_arg = [{"Dim": er.dim, "N_points": er.n_points,
               "Epsilon": er.epsilon,
               "Point_Generation_Method": er.point_generation_method,
               "Spanner_Method": er.spanner_method, "Spanner_Edges": er.spanner_edges,
               "Sparseness": er.sparseness} for r in results for er in r]

    df = pd.DataFrame(df_arg)
    df.to_pickle("blind_random_spanner_results_pandas-quasi.pkl")

    for r in results:
        for er in r:
            print(er)
