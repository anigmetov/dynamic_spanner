#!/usr/bin/env python3

import sys
import numpy as np
from utils_euclidean import *

def save_to_files(point_method, dim, n_points, point_file_name, matrix_file_name = None):
    if matrix_file_name == None:
        matrix_file_name = "dist_matrix_" + point_file_name
    pm_dic = {"unif": get_uniform_points,
              "normal": get_points,
              "clustered": get_clustered_points,
              "exp": get_exponential_points}
    if not point_method in pm_dic:
        print("Unknown point_method")
        exit(1)
    points = pm_dic[point_method](n_points, dim)
    np.savetxt(point_file_name, points, "%15.7f")
    dist_matr = distance_matrix(points)
    np.savetxt(matrix_file_name, dist_matr, "%15.7f")


if __name__ == "__main__":


    for n_points in range(50, 850, 50):
        for dim in [2, 3, 10]:
            for point_method in ["unif", "normal", "exp", "clustered"]:
                point_file_name = "points_{}_{}_{}.txt".format(dim, n_points, point_method)
                np.random.seed(1)
                save_to_files(point_method, dim, n_points, point_file_name)
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


