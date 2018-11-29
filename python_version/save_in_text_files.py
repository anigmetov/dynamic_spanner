#!/usr/bin/env python3

import sys
import numpy as np
from utils_euclidean import *
import joblib as jl


def save_to_files(point_method, dim, n_points, point_file_name=None, matrix_file_name=None, norm = None):
    if point_file_name == None:
        point_file_name = "points_{}_{}_{}.txt".format(dim, n_points, point_method)
    if matrix_file_name == None:
        matrix_file_name = "dist_matrix_" + point_file_name
    pm_dic = {"unif": get_uniform_points,
              "normal": get_points,
              "clustered": get_clustered_points,
              "exp": get_exponential_points,
              "beta_05_05": get_beta_05_05,
              "beta_5_1": get_beta_5_1,
              "beta_1_3": get_beta_1_3,
              "beta_2_2": get_beta_2_2,
              "beta_2_5": get_beta_2_5,
              }
    if not point_method in pm_dic:
        print("Unknown point_method")
        exit(1)
    points = pm_dic[point_method](n_points, dim)
    np.savetxt(point_file_name, points, "%15.7f")
    dist_matr = distance_matrix(points)
    np.savetxt(matrix_file_name, dist_matr, "%15.7f")

def save_to_files(point_method, dim, n_points, norm = 2.0, point_file_name=None, matrix_file_name=None):
    if point_file_name == None:
        point_file_name = "points_{}_{}_{}.txt".format(dim, n_points, point_method)
    if matrix_file_name == None:
        matrix_file_name = "dist_matrix_norm_{}_{}".format(norm,  point_file_name)
    pm_dic = {"unif": get_uniform_points,
              "normal": get_points,
              "clustered": get_clustered_points,
              "exp": get_exponential_points,
              "beta_05_05": get_beta_05_05,
              "beta_5_1": get_beta_5_1,
              "beta_1_3": get_beta_1_3,
              "beta_2_2": get_beta_2_2,
              "beta_2_5": get_beta_2_5,
              }
    if not point_method in pm_dic:
        print("Unknown point_method")
        exit(1)
    points = pm_dic[point_method](n_points, dim)
    np.savetxt(point_file_name, points, "%15.7f")
    dist_matr = distance_matrix(points, norm)
    np.savetxt(matrix_file_name, dist_matr, "%15.7f")

def save_unif_files(dim, n_points, max_coord):
    point_file_name = "points_uniform_dim_{}_max_coord_{}_npoints_{}.txt".format(dim, max_coord, n_points)
    matrix_file_name = "dist_matrix_" + point_file_name
    points = get_uniform_points(n_points, dim, max_coord)
    np.savetxt(point_file_name, points, "%13.7f")
    dist_matr = distance_matrix(points)
    np.savetxt(matrix_file_name, dist_matr, "%13.7f")



if __name__ == "__main__":


    if False:
        results = jl.Parallel(n_jobs=-1)(
            jl.delayed(save_unif_files)(dim, n_points, max_coord)
            for dim in [2, 3, 4]
            # for n_points in range(50, 500, 50)
            for n_points in [ 500 ]
            for max_coord in [10.0 ** i for i in range(0, 4)])


    if True:
        results = jl.Parallel(n_jobs=-1)(
            jl.delayed(save_to_files)(point_method, dim, n_points, norm)
            for dim in [2, 3, 4, 5]
            for norm in [1.0, 3.0, np.inf]
            for n_points in range(50, 550, 50)
            for point_method in ["unif", "normal"]
            # for point_method in ["unif", "normal", "exp", "clustered", "beta_05_05",
            #  "beta_5_1", "beta_1_3", "beta_2_2", "beta_2_5"]
        )

        # for n_points in range(50, 850, 50):
        #     for dim in [4, 5, 6]:
        #         # for point_method in ["beta_05_05", "beta_5_1", "beta_1_3", "beta_2_2", "beta_2_5"]:
        #         # for point_method in ["unif", "normal", "exp", "clustered"]:
        #         for point_method in ["unif", "normal", "exp", "clustered", "beta_05_05", "beta_5_1", "beta_1_3", "beta_2_2", "beta_2_5"]:
        #             np.random.seed(1)
        #             save_to_files(point_method, dim, n_points, point_file_name)
        sys.exit(0)

    if len(sys.argv) < 5:
        print("Parameters: point_method(unif, normal, exp, clustered) dim n_points points_file [matrix_file]")
        exit(1)
    arg_idx = 1
    point_method = sys.argv[arg_idx]
    arg_idx += 1
    dim = int(sys.argv[arg_idx])
    arg_idx += 1
    n_points = int(sys.argv[arg_idx])
    arg_idx += 1
    point_file_name = sys.argv[arg_idx]
    arg_idx += 1
    if len(sys.argv) == arg_idx:
        matrix_file_name = ""
    else:
        matrix_file_name = sys.argv[arg_idx]

    save_to_files(point_method, dim, n_points, point_file_name, matrix_file_name)
