#!/usr/bin/env python3

import os
import os.path
import sys

import numpy as np
import dionysus as d
import datetime as dt


class MyDiagram:
    def __init__(self, dgm_data, id, label=None):
        self.dgm_data = dgm_data
        self.id = id
        self.label = label



def read_mcgill_diagrams(diag_dir_root, dim):
    diag_dirs = [os.path.join(diag_dir_root, f) for f in os.listdir(diag_dir_root)
                 if f.endswith('Ply')]
    dim_ext = '.%d.txt' % dim
    dgms = []
    dgm_id = 0
    for diag_dir in diag_dirs:
        diag_files = [os.path.join(diag_dir, f) for f in os.listdir(diag_dir)
                      if f.endswith(dim_ext)]
        for dgm_f in diag_files:
            pts = []
            with open(dgm_f, 'r') as dgm_file:
                for line in dgm_file:
                    x, y = [float(x) for x in line.split()]
                    pts.append((x, y))
            dgms.append(MyDiagram(d.Diagram(pts[:]), dgm_id, dgm_f))
            dgm_id += 1
            # print("Read file %s, number of points %d" % (dgm_f, len(pts)))
    print("Read %d diagrams" % len(dgms))
    return dgms


def get_distance_tuple(i, j, dgm1, dgm2, q, delta):
    # print("i = %d, j = %d, q = %d, fname1 = %s, fname2 = %s, started computation" % (i, j, q, dgm1.label, dgm2.label))
    start_time = dt.datetime.now()
    res = d.wasserstein_distance(dgm1.dgm_data, dgm2.dgm_data, q=q, delta=delta)
    end_time = dt.datetime.now()
    elapsed_time = (end_time - start_time).microseconds
    # print("i = %d, j = %d, q = %d, fname1 = %s, fname2 = %s, d = %f, time = %d microsec, finished computation" % (i, j, q, dgm1.label, dgm2.label, res, elapsed_time))
    return (i, j, elapsed_time, elapsed_time)
    

if __name__ == "__main__":
    # if len(sys.argv) < 2:
    #     print("dim\n")
    #     exit(0)
    # for dim in [2, 1, 0]:
    for dim in [2,0]:
    # dim = int(sys.argv[1])
        for delta in [0.01, 0.5, 0.2, 0.1]:
        # for delta in [0.5, 0.2, 0.1, 0.01]:
            diag_dir_root = '/home/anigmetov/Downloads/mcgill_shape_benchmark_dgms_no_duplicates/'
            dgms = read_mcgill_diagrams(diag_dir_root, dim)
            n_pts = len(dgms)
            print("n_pts = %d\n" % n_pts)
            total_distances = float(n_pts * (n_pts - 1) / 2)
            data = np.array(dgms[:])
            times = {}
            for q in [1,2,3]:
                for attempt in range(5):
            # for q in [3]:

                    times[q] = 0
                    dist_matr_fname = "dist_timing_matrix_mcgill_original_q_{}_dim_{}_delta_{}_attempt_{}.txt".format(q, dim, delta, attempt)
                    dist_matrix = np.zeros((n_pts, n_pts))
                    for i in range(n_pts):
                        for j in range(i+1, n_pts):
                            (i, j, wd, et) = get_distance_tuple(i, j, dgms[i], dgms[j], q, delta)
                            dist_matrix[i][j] = wd
                            dist_matrix[j][i] = wd
                            times[q] += et
                    print("ELAPSED TIME = {} ms = {} sec for q = {}, dim = {}, delta = {}".format(times[q], times[q] / 1e6, q, dim, delta))

                    with open(dist_matr_fname, 'w') as dmf:
                        dmf.write("%d\n" % n_pts)
                        for i in range(n_pts):
                            for j in range(0, n_pts):
                                dmf.write("%f\n" % dist_matrix[i][j])
                        dmf.write("\n%f seconds\n" % (times[q] / 1e6))
