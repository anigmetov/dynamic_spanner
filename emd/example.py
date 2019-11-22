#!/usr/bin/env python

from PIL import Image
import numpy as np
from emd import emd

def get_pic_support(pic_shape):
    result = []
    for i in range(pic_shape[0]):
        for j in range(pic_shape[1]):
            result.append([float(i),float(j)])
    result = np.asarray(result, dtype=np.float64)
    return result

def read_image(fname):
    pic = Image.open(fname)
    pic_weights = np.asarray(pic.convert("L"), dtype=np.float64)
    pic_weights = pic_weights / np.sum(pic_weights)
    pic_shape = pic_weights.shape
    n_pixels = pic_shape[0] * pic_shape[1]
    pic_weights = pic_weights.reshape((n_pixels, 1))
    pic_support = get_pic_support(pic_shape)
    return (pic_weights, pic_support)


if __name__ == "__main__":
    a_weights, A = read_image("00001.ppm")
    b_weights, B = read_image("00002.ppm")
    print("read images")
    res = emd(A, B, a_weights, b_weights)
    print("distance = ", res)
