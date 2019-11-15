#!/usr/bin/env python

import numpy as np

if __name__ == "__main__":
    flist = [ "matrix-try" ]

    n_procs = 10
    n_categories = 10
    n_images_per_cat = 20

    n_images = n_categories * n_images_per_cat
    n_images = 13
    n_pairs = n_images * (n_images - 1) / 2

    d_matr = np.zeros((n_pairs, n_pairs), dtype = np.float32)
    t_matr = np.zeros((n_pairs, n_pairs), dtype = np.float32)

    for fname in flist:
        with open(fname, 'r') as f:
            for line in f:
                a = line.split('\t')
                i = int(a[0])
                j = int(a[1])
                d_matr[i][j] = d_matr[j][i] = float(a[4])
                t_matr[i][j] = t_matr[j][i] = float(a[5])

    np.savetxt("matrix-try-dist", d_matr)
    np.savetxt("matrix-try-time", t_matr)
