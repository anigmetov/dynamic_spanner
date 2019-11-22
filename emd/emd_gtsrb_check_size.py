#!/usr/bin/env python

from timeit import default_timer as timer
import sys

import glob
import random


from PIL import Image
import numpy as np
from emd import emd

from emd_utils import *

if __name__ == "__main__":
    random.seed(42)
    n_images = 300
    im_files = read_gtsrb_file_list(n_images)
    a = []
    for fname in im_files:
        r = load_picture(fname)
        a.append(r)

    a.sort()

    for t in a:
        print(t)

