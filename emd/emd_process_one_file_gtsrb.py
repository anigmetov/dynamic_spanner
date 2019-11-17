#!/usr/bin/env python

from timeit import default_timer as timer
import sys

from PIL import Image
import numpy as np
from emd import emd

# return flattened list of pixel coordinates (support of measure)
def pic_positions(shape):
    result = []
    for i in range(shape[0]):
        for j in range(shape[1]):
            result.append([float(i),float(j)])
    result = np.asarray(result, dtype=np.float64)
    return result

def load_picture(fname):
    a = Image.open(fname)
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


if __name__ == "__main__":
    fname_in = sys.argv[1]
    part_idx = fname_in.split(".")[-1]
    fname_out = sys.argv[2] + "." + part_idx
    compute_from_file(fname_in, fname_out)
