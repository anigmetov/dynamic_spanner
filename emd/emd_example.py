#!/usr/bin/env python

from PIL import Image
import numpy as np
from emd import emd

def weights(pic_size):
    result = []
    for i in range(pic_size):
        for j in range(pic_size):
            result.append([float(i),float(j)])
    result = np.asarray(result, dtype=np.float64)
    return result

pic_size = 32
n_pixels = pic_size * pic_size

a = Image.open("bismarck.jpg")
a = a.resize((pic_size,pic_size))
a_weights = np.asarray(a.convert("L"), dtype=np.float64)
a_weights = a_weights / np.sum(a_weights)

b = Image.open("kipling.jpg")
b = b.resize((pic_size,pic_size))
b_weights = np.asarray(b.convert("L"), dtype=np.float64)
b_weights = b_weights / np.sum(b_weights)

a_weights = a_weights.reshape((n_pixels, 1))
b_weights = b_weights.reshape((n_pixels, 1))

print("read images")

A = weights(pic_size)
B = weights(pic_size)

# print(A)
# print(B)
# print(a_weights)
# print(b_weights)

res = emd(A, B, a_weights, b_weights)
print("distance = ", res)
