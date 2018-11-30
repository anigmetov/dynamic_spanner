#!/usr/bin/env python3

import sys
import numpy as np

def get_lp_distance(a, b, norm = 2.0):
    return np.linalg.norm(a - b, ord=norm)

def get_points(n_points, dim):
    return np.random.randn(n_points, dim)


def get_normal_points(n_points, dim, scale = 1.0):
    return np.random.normal(scale = scale, size = (n_points, dim))

def get_uniform_points(n_points, dim, max_coord=100.0):
    return np.random.uniform(0.0, max_coord, (n_points, dim))


def get_exponential_points(n_points, dim):
    max_power = 25  # max(n_points / 2.0, 10.0)
    # print(f"max_power = {max_power}")
    shape = (n_points, dim)
    res = np.random.uniform(1.0, max_power, shape)
    return np.power(2.0, res)


def get_clustered_points(n_points, dim, n_clusters = -1):
    if n_clusters == -1:
        n_clusters = max(5, int(n_points / 50))
    assert(n_points % n_clusters == 0)
    points_per_cluster = int(n_points / n_clusters)
    a = get_uniform_points(n_clusters, dim, 10000.0)
    cluster_centres = np.repeat(a, points_per_cluster, axis = 0)
    # print("********************************************************************************")
    # print(a)
    # print("********************************************************************************")
    # print(cluster_centres)
    # print("********************************************************************************")
    diffs = np.random.normal(scale = 10.0, size = (n_points, dim))
    return diffs + cluster_centres


def get_beta_05_05(n_points, dim):
    return np.random.beta(0.5, 0.5, (n_points, dim))

def get_beta_5_1(n_points, dim):
    return np.random.beta(5, 1, (n_points, dim))

def get_beta_1_3(n_points, dim):
    return np.random.beta(1, 3, (n_points, dim))

def get_beta_2_2(n_points, dim):
    return np.random.beta(2, 2, (n_points, dim))

def get_beta_2_5(n_points, dim):
    return np.random.beta(2, 5, (n_points, dim))

def distance_matrix(points, norm = 2.0):
    n_points = len(points)
    result = np.zeros((n_points, n_points))
    for i in range(n_points):
        for j in range(n_points):
            if i <= j:
                continue
            result[i][j] = result[j][i] = get_lp_distance(points[i, :], points[j, :], norm)
    return result


def distances_list(points):
    n_points = len(points)
    dm = distance_matrix(points)
    return [(dm[i][j], i, j) for i in range(n_points) for j in range(i + 1, n_points)]


def sorted_distances(points):
    return sorted(distances_list(points))

