#!/usr/bin/env python3

class ExperimentResult:
    def __init__(self, dim, n_points, eps, points_gen_method, spanner_method, spanner_edges):
        self.dim = dim
        self.n_points = n_points
        self.epsilon = eps
        self.total_n_pairs = n_points * (n_points - 1) / 2.0
        self.spanner_method = spanner_method
        self.spanner_edges = spanner_edges
        self.sparseness = spanner_edges / self.total_n_pairs
        self.point_generation_method = points_gen_method

    def __str__(self):
        return "{};{};{};{};{};{};{}".format(self.spanner_method, self.dim, self.n_points, self.epsilon,
                                             self.point_generation_method, self.spanner_edges, self.sparseness)

