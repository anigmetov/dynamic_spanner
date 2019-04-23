#!/usr/bin/env python3

import os
import os.path
import sys
import re

import numpy as np

def get_time_from_last_line(fname):
    with open(dist_matr_fname, 'r') as dmf:
        last_line = ""
        for s in dmf:
            last_line = s
        reg_res = re.match("(.*) seconds")
        if reg_res:
            return float(reg_res.group(1))
        else:
            return 0.0


if __name__ == "__main__":
    for dim in [2,0]:
        for delta in [0.01, 0.5, 0.2, 0.1]:
            for q in [1,2,3]:
                total_time = 0.0
                n_attempts = 5
                for attempt in range(n_attempts):
                    dist_matr_fname = "dist_timing_matrix_mcgill_original_q_{}_dim_{}_delta_{}_attempt_{}.txt".format(q, dim, delta, attempt)
                    time_in_file = get_time_from_last_line(fname)
                    total_time += time_in_file
                avg_time = total_time / n_attempts
