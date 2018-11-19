#!/usr/bin/env python3

import sys
import numpy as np
from utils_euclidean import *

if __name__ == "__main__":
    if len(sys.argv) < 5:
        print("Parameters: point_method(unif or normal) dim n_points points_file [matrix_file]")
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
        matrix_file_name = "dist_matrix_" + point_file_name
    else:
        matrix_file_name = sys.argv[arg_idx]

    pm_dic = { "unif" : get_uniform_points, "normal" : get_points }
    if not point_method in pm_dic:
        print("Unknown point_method")
        exit(1)
    points = pm_dic[point_method](n_points, dim)
    np.savetxt(point_file_name, points, "%15.7f")
    dist_matr = distance_matrix(points)
    np.savetxt(matrix_file_name, dist_matr, "%15.7f")