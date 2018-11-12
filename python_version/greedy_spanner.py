#!/usr/bin/env python3

import networkx as nx
import sys
import matplotlib.pyplot as plt
import joblib as jl

from spanner_experiment_result import ExperimentResult
from weakly_sorted import *

from utils_euclidean import *

import pandas as pd

inf_dist = 10000000000


def upper_bound(g, i, j):
    try:
        return nx.shortest_path_length(g, i, j, weight='weight')
    except nx.NetworkXNoPath:
        return inf_dist


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


def run_experiment(dim, n_points, epsilon, points_generator, ps_gen_args):
    points = points_generator(n_points, dim, *ps_gen_args)
    dist_list = distances_list(points)

    result = []

    gs = greedy_spanner(points, epsilon, sorted(dist_list))
    n_gs_edges = gs.number_of_edges()

    result.append(ExperimentResult(dim, n_points, epsilon, points_generator.__name__, "greedy", n_gs_edges))

    print(result[-1])

    random_dist_list = dist_list[:]
    random.shuffle(random_dist_list)
    gs_random = greedy_spanner(points, epsilon, random_dist_list)

    n_gs_random_edges = gs_random.number_of_edges()
    result.append(ExperimentResult(dim, n_points, epsilon, points_generator.__name__, "random", n_gs_random_edges))
    print(result[-1])

    sd = quasi_sorted_distances(points, dist_list)
    gs_quasi = greedy_spanner(points, epsilon, sd.qsorted_distances)
    n_gs_quasi_edges = gs_quasi.number_of_edges()
    result.append(ExperimentResult(dim, n_points, epsilon, points_generator.__name__, "quasi", n_gs_quasi_edges))
    print(result[-1])

    sys.stdout.flush()

    return result


if __name__ == "__main__":
    # min_dist = 0.01
    # max_dist = 3.0
    # print(get_buckets(min_dist, max_dist))

    # n_pointses = [100]
    # epsilons = [0.1]
    # dims = [2]
    # ps_gen_methods = [get_points, get_exponential_points, get_uniform_points]


    n_pointses = [800]
    dims = [2, 3, 4]
    ps_gen_methods = [get_points, get_uniform_points, get_exponential_points]
    ps_gen_args = [[], [10.0], []]
    epsilons = [0.01, 0.1, 0.2, 0.5, 2.0]

    all_results = jl.Parallel(n_jobs=5)(
        jl.delayed(run_experiment)(dim, n_points, epsilon, ps_gen_method, ps_arg)
        for dim in dims
        for n_points in n_pointses
        for epsilon in epsilons
        for (ps_gen_method, ps_arg) in zip(ps_gen_methods, ps_gen_args))

    print("################################################################################")

    df_arg = [{"Dim": er.dim, "N_points": er.n_points,
               "Epsilon": er.epsilon,
               "Point_Generation_Method": er.point_generation_method,
               "Spanner_Method": er.spanner_method, "Spanner_Edges": er.spanner_edges,
               "Sparseness": er.sparseness} for r in all_results for er in r]

    df = pd.DataFrame(df_arg)
    df.to_pickle("greedy_spanner_results_pandas-800.pkl")

    for r in all_results:
        for er in r:
            print(er)

    sys.exit(0)
