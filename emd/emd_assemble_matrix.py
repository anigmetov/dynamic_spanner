#!/usr/bin/env python

import numpy as np

from emd_utils import *

if __name__ == "__main__":
    flist = [ "matrix-try" ]

    n_procs = 10
    n_categories = 10
    n_images_per_cat = 20

    im_size = 35

    flist = [ "pair_file_nc_%d_i_%d_np_%d.%d.matrix.%d.%d" % (n_categories, n_images_per_cat, n_procs, i, im_size, i)
            for i in range(n_procs) ]

    n_images = n_categories * n_images_per_cat
    n_pairs = n_images * (n_images - 1) / 2

    d_matr = np.zeros((n_images, n_images), dtype = np.float32)
    t_matr = np.zeros((n_images, n_images), dtype = np.float32)

    for fname in flist:
        print(fname)
        with open(fname, 'r') as f:
            for line in f:
                a = line.split('\t')
                i = int(a[0])
                j = int(a[1])
                d_matr[i][j] = d_matr[j][i] = float(a[4])
                t_matr[i][j] = t_matr[j][i] = float(a[5])

    d_matr_fname = "dist_matr_nc_%d_i_%d_size_%d.txt" % (n_categories, n_images_per_cat, im_size)
    t_matr_fname = "time_matr_nc_%d_i_%d_size_%d.txt" % (n_categories, n_images_per_cat, im_size)

    save_matrix(d_matr_fname, d_matr)
    save_matrix(t_matr_fname, t_matr)

    # # np.savetxt("matrix-try-dist", d_matr)
    # np.savetxt("matrix-try-time", t_matr)
