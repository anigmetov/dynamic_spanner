#!/usr/bin/env python3

import os
import os.path
import sys

import numpy as np
import dionysus as d

wasser_dist_calls = 0
wasser_dist_dict = {}
q = 1


class MyDiagram:
    def __init__(self, dgm_data, id, label=None):
        self.dgm_data = dgm_data
        self.id = id
        self.label = label


def wasser_dist(dgm1, dgm2, q):
    global wasser_dist_calls
    global wasser_dist_dict
    if dgm1.id == dgm2.id:
        return 0.0
    key1 = (dgm1.id, dgm2.id)
    key2 = (dgm2.id, dgm1.id)
    if key1 in wasser_dist_dict:
        return wasser_dist_dict[key1]
    else:
        wasser_dist_calls += 1
        # print("wasser_dist called, %d, %d, calling Hera" % (dgm1.id, dgm2.id))
        result = d.wasserstein_distance(dgm1.dgm_data, dgm2.dgm_data, q=q,
                                        delta=0.01)
        # print("Hera returned")
        wasser_dist_dict[key1] = result
        wasser_dist_dict[key2] = result
        return result


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
            print("Read file %s, number of points %d" % (dgm_f, len(pts)))
    print("Read %d diagrams" % len(dgms))
    return dgms


def get_distance_tuple(i, j, dgm1, dgm2, q):
    print("i = %d, j = %d, q = %d, fname1 = %s, fname2 = %s, started computation" % (i, j, q, dgm1.label, dgm2.label))
    res = wasser_dist(dgm1, dgm2, q)
    print("i = %d, j = %d, q = %d, fname1 = %s, fname2 = %s, d = %f, finished computation" % (i, j, q, dgm1.label, dgm2.label, res))
    return (i, j, res)
    

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("dim\n")
        exit(0)
    dim = int(sys.argv[1])
    diag_dir_root = '/home/anigmetov/Downloads/mcgill_shape_benchmark_dgms_no_duplicates/'
    dgms = read_mcgill_diagrams(diag_dir_root, dim)
    n_pts = len(dgms)
    print("n_pts = %d\n" % n_pts)
    total_distances = float(n_pts * (n_pts - 1) / 2)
    data = np.array(dgms[:])
    for q in [1, 2, 3]:
        dist_matr_fname = "dist_matrix_mcgill_original_no_parallel_q_{}_dim_{}.txt".format(q, dim)
        dist_matr_fname_npy = "dist_matrix_mcgill_original_no_parallel_q_{}_dim_{}.npy".format(q, dim)
        wasser_dist_calls = 0
        wasser_dist_dict = {}
        dist_matrix = np.zeros((n_pts, n_pts))
        dist_list = [ get_distance_tuple(i, j, dgms[i], dgms[j], q) for i in range(n_pts) for j in range(i+1, n_pts)] 
        for t in dist_list:
            i, j, wd = t
            dist_matrix[i][j] = wd
            dist_matrix[j][i] = wd

        np.save(dist_matr_fname_npy, dist_matrix)

        with open(dist_matr_fname, 'w') as dmf:
            dmf.write("%d\n" % n_pts)
            for i in range(n_pts):
                for j in range(n_pts):
                    dmf.write("%f\n" % dist_matrix[i][j])
