#!/usr/bin/env python

from timeit import default_timer as timer
import sys

import glob
import random


from PIL import Image
import numpy as np
from emd import emd


def read_gtsrb_file_list(n_images):
    im_files_all = glob.glob("traffic_signs/*.ppm")
    im_files = random.sample(im_files_all, n_images)
    return im_files


def resize_pic(pic, coeff):
    new_width = pic.width / coeff
    new_height = pic.height / coeff
    return pic.resize((new_width, new_height), Image.BICUBIC)



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


def save_matrix(fname, a):
    n_rows = a.shape[0]
    n_cols = a.shape[1]
    n_entries = n_rows * n_cols
    with open(fname, 'w') as f:
        f.write("%d\n" % n_rows)
        for i in range(n_rows):
            for j in range(n_cols):
                f.write("%f\n" % a[i][j])

# return flattened list of pixel coordinates (support of measure)
def pic_positions(shape):
    result = []
    for i in range(shape[0]):
        for j in range(shape[1]):
            result.append([float(i),float(j)])
    result = np.asarray(result, dtype=np.float64)
    return result


def resize_pic(pic, coeff):
    new_width = pic.width / coeff
    new_height = pic.height / coeff
    return pic.resize((new_width, new_height), Image.BICUBIC)


def load_picture(fname):
    a = Image.open(fname)

    a = resize_pic(a, 2)

    # get grayvalues
    a_weights = np.asarray(a.convert("L"), dtype=np.float64)
    # normalize
    a_weights = a_weights / np.sum(a_weights)

    # get image support
    a_positions = pic_positions(a_weights.shape)

    # flatten
    n_pixels = a_weights.shape[0] * a_weights.shape[1]
    a_weights = a_weights.reshape((n_pixels, 1))

    return (a_weights, a_positions)


def compute_distance(fname_a, fname_b):
    weights_a, positions_a = load_picture(fname_a)
    weights_b, positions_b = load_picture(fname_b)
    start = timer()
    dval = emd(positions_a, positions_b, weights_a, weights_b)
    end = timer()
    elapsed = end - start
    return (dval, elapsed)


def compute_from_file(fname_in, fname_out):
    with open(fname_in, 'r') as pair_file, open(fname_out, 'w') as matrix_file:
        for line in pair_file:
            fname_a, fname_b, i, j = line.split(' ')
            i = int(i)
            j = int(j)
            dval, elapsed = compute_distance(fname_a, fname_b)
            matrix_file.write("%d\t%d\t%s\t%s\t%f\t%f\n" % (i, j, fname_a, fname_b, dval, elapsed))
            matrix_file.flush()
