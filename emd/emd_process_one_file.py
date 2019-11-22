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


def compute_distance(fname_a, fname_b, pic_size):
    pic_shape = (pic_size, pic_size)
    positions = pic_positions(pic_shape)
    weights_a = load_picture(fname_a, pic_size)
    weights_b = load_picture(fname_b, pic_size)
    start = timer()
    dval = emd(positions, positions, weights_a, weights_b)
    end = timer()
    elapsed = end - start
    return (dval, elapsed)


def compute_from_file(fname_in, fname_out, pic_size):
    with open(fname_in, 'r') as pair_file, open(fname_out, 'w') as matrix_file:
        for line in pair_file:
            fname_a, fname_b, i, j = line.split(' ')
            i = int(i)
            j = int(j)
            dval, elapsed = compute_distance(fname_a, fname_b, pic_size)
            matrix_file.write("%d\t%d\t%s\t%s\t%f\t%f\n" % (i, j, fname_a, fname_b, dval, elapsed))


if __name__ == "__main__":
    fname_in = sys.argv[1]
    part_idx = fname_in.split(".")[-1]
    fname_out = sys.argv[2] + "." + part_idx
    pic_size = int(sys.argv[3])
    compute_from_file(fname_in, fname_out, pic_size)
