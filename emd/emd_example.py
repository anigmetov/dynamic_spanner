#!/usr/bin/env python

import glob
from timeit import default_timer as timer
import random

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
    print("positions =  ", result.shape)
    return result

def load_picture(fname, pic_size = None):
    a = Image.open(fname)
    if pic_size:
        a = a.resize((pic_size,pic_size))
    # get grayvalues
    a_weights = np.asarray(a.convert("L"), dtype=np.float64)
    # normalize
    a_weights = a_weights / np.sum(a_weights)
    # flatten
    n_pixels = a_weights.shape[0] * a_weights.shape[1]
    a_weights = a_weights.reshape((n_pixels, 1))

    return a_weights

def prepare_caltech(n_categories, n_images_per_cat):
    cat_dirs_all = glob.glob("caltech_101/*")
    cat_dirs = random.sample(cat_dirs_all, n_categories)
    im_files = []
    for cat_dir in cat_dirs:
        im_files_all = glob.glob(cat_dir + "/*.jpg")
        im_files.extend(random.sample(im_files_all, n_images_per_cat))
    return im_files





if __name__ == "__main__":

    random.seed(1)

    n_categories = 12
    n_images_per_cat = 30

    im_files = prepare_caltech(n_categories, n_images_per_cat)

    elapsed_total = 0.0

    pic_size = 40
    pic_shape = (pic_size, pic_size)
    positions = pic_positions(pic_shape)

    out_fname = "caltech_cats_%d_pic_per_cat_%d_size_%d.txt" % (n_categories,
            n_images_per_cat, pic_size)

    log_fname = "log_caltech_cats_%d_pic_per_cat_%d_size_%d.txt" % (n_categories,
            n_images_per_cat, pic_size)

    with open(out_fname, 'w') as f, open(log_fname, 'w') as lf: 
        for i in range(len(im_files)):
            weights_i = load_picture(im_files[i], pic_size)
            for j in range(0, i):
                weights_j = load_picture(im_files[j], pic_size)
                start = timer()
                dval = emd(positions, positions, weights_i, weights_j)
                end = timer()
                elapsed_ij = end - start
                elapsed_total += elapsed_ij
                f.write(("%d %d %f\n" %(i, j, dval)))
                lf.write("%d %d %s %s dist = %f elapsed = %f sec\n" % (i, j, im_files[i], im_files[j], dval, elapsed_ij))
        lf.write("total elapsed = %f sec\n" % elapsed_total)
