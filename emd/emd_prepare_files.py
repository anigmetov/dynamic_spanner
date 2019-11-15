#!/usr/bin/env python

import glob
import random

def read_caltech_file_list(n_categories, n_images_per_cat):
    cat_dirs_all = glob.glob("caltech_101/*")
    cat_dirs = random.sample(cat_dirs_all, n_categories)
    im_files = []
    for cat_dir in cat_dirs:
        im_files_all = glob.glob(cat_dir + "/*.jpg")
        im_files.extend(random.sample(im_files_all, n_images_per_cat))
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

    random.seed(1)

    n_categories = 2
    n_images_per_cat = 4

    n_procs = 4

    im_files = read_caltech_file_list(n_categories, n_images_per_cat)
    prepare_file_lists_for_parallel(im_files, n_procs, 'pair_file')


