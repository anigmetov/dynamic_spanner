#!/usr/bin/env python

from PIL import Image
import numpy as np

a = Image.open("bismarck.jpg")
a = a.resize(128,128)
a_weights = np.asarray(a.convert("L"), dtype=np.float32)
a_weights = a_weights / np.sum(a_weights)


