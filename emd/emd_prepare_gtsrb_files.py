#!/usr/bin/env python

import glob
import random

def read_gtsrb_file_list(n_images):
    im_files_all = glob.glob("traffic_signs/*.ppm")
    im_files = random.sample(im_files_all, n_images)
    return im_files

def split(a, n):
    k, m = divmod(len(a), n)
    return (a[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in xrange(n))


def prepare_file_lists_for_parallel(im_files, n_procs, file_start):
    all_pairs = [ (i, j) for i in range(len(im_files)) for j in range(0, i) ]
    all_pairs_split = split(all_pairs, n_procs)
    for i, part in enumerate(all_pairs_split):
        fname = file_start + "." + str(i)
        with open(fname, 'w') as f:
            for x, y in part:
                f.write("%s %s %d %d\n" % (im_files[x], im_files[y], x, y))


if __name__ == "__main__":

    random.seed(42)

    n_images = 400

    n_procs = 10

    fname_out = "gtsrb_pair_file_i_%d_np_%d" % (n_images, n_procs)

    im_files = read_gtsrb_file_list(n_images)
    prepare_file_lists_for_parallel(im_files, n_procs, fname_out)
