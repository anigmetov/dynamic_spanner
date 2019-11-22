#!/usr/bin/env python

import random
import numpy as np

from emd_utils import *

if __name__ == "__main__":
    flist = [ "matrix-try" ]

    n_procs = 10
    n_images = 300
    n_pairs = n_images * (n_images - 1) / 2

    flist = [ "gtsrb_pair_file_i_%d_np_%d.%d.matrix.%d" % (n_images, n_procs, i, i) for i in range(n_procs) ]
    flist_1 = [ "gtsrb_missing_.%d.matrix.%d" % (i, i) for i in range(n_procs) ]
    flist.extend(flist_1)


    d_matr = np.zeros((n_images, n_images), dtype = np.float32)
    t_matr = np.zeros((n_images, n_images), dtype = np.float32)

    image_fnames = [ "" for i in range(n_images) ]

    for fname in flist:
        print(fname)
        with open(fname, 'r') as f:
            for line in f:
                a = line.split('\t')
                i = int(a[0])
                j = int(a[1])
                image_fnames[i] = a[2]
                image_fnames[j] = a[3]
                d_matr[i][j] = d_matr[j][i] = float(a[4])
                t_matr[i][j] = t_matr[j][i] = 1000000 * float(a[5])

    computed_distances = np.sum(d_matr > 0.0, 1) + 1

    complete_indices = np.where(computed_distances == n_images)
    tot_incomplete_indices = np.where(computed_distances == 1)
    partial_indices = np.where((computed_distances > 1) & (computed_distances < n_images))

    print("complete: ", complete_indices[0].size)

    partial_image_fnames = [ (image_fnames[i], i) for i in partial_indices[0] ]

    # # print("partial: ", partial_image_fnames)

    small_partial_image_fnames = []
    max_size = 900
    for pif, pif_idx in partial_image_fnames:
        w, p = load_picture(pif)
        if p.shape[0] > max_size:
            print(pif, " is too large: ", p.shape)
        else:
            small_partial_image_fnames.append((pif, pif_idx))

    small_partial_indices = [ pif_idx for (pif, pif_idx) in small_partial_image_fnames ]
    aaa = [ x for x in complete_indices[0] ]
    aaa.extend(small_partial_indices)

    aaa = random.sample(aaa, 250)

    aaa.sort()
    bbb = set(aaa)

    print(aaa)
    print(len(aaa), len(bbb))


    final_indices = np.asarray(aaa)


    
    # print("small partial: ", small_partial_image_fnames)

    # missing_pairs = []
    # for pif, pif_idx in small_partial_image_fnames:
    #     for j in range(0, pif_idx):
    #         if d_matr[pif_idx][j] <= 0.0:
    #             missing_pairs.append((pif_idx, j))

    # print("missing pairs: ", len(missing_pairs))

    # add_pair_prefix = "gtsrb_missing_"

    # n_procs = 10
    # all_pairs_split = split(missing_pairs, n_procs)
    # for i, part in enumerate(all_pairs_split):
    #     fname = add_pair_prefix + "." + str(i)
    #     with open(fname, 'w') as f:
    #         for x, y in part:
    #             f.write("%s %s %d %d\n" % (image_fnames[x], image_fnames[y], x, y))


    # for i, j in missing_pairs:
    #     print("%s %s %d %d" % (image_fnames[i], image_fnames[j], i, j))


    compl_ind = np.ix_(final_indices, final_indices)
    d_matr = d_matr[compl_ind]
    t_matr = t_matr[compl_ind]

    total_time = np.sum(t_matr, (0,1)) * 0.5
    print("total time", total_time)

    n_images = final_indices.size

    d_matr_fname = "dist_matr_images_%d.txt" % (n_images)
    t_matr_fname = "time_matr_images_%d.txt" % (n_images)

    save_matrix(d_matr_fname, d_matr)
    save_matrix(t_matr_fname, t_matr)

    # # np.savetxt("matrix-try-dist", d_matr)
    # np.savetxt("matrix-try-time", t_matr)
