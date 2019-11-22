#!/usr/bin/env python

import glob
import random

from emd_utils import *

if __name__ == "__main__":

    random.seed(42)

    n_images = 300

    n_procs = 10

    fname_out = "gtsrb_pair_file_i_%d_np_%d" % (n_images, n_procs)

    im_files = read_gtsrb_file_list(n_images)
    prepare_file_lists_for_parallel(im_files, n_procs, fname_out)
