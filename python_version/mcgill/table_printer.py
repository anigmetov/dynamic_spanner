#!/usr/bin/env python3

import os
import os.path
import glob
import sys
import re
from math import floor


def sum_file(fname):
    result = 0
    with open(fname, 'r') as f:
        first_line = True
        for s in f:
            if first_line:
                first_line = False
                continue
            result += float(s)
    return floor(result / 2e6)

if __name__ == "__main__":
    epss = ["0.1", "0.2", "0.5"]
    qs = [1,2,3]
    dims = [0]

    final_table = {}

    for dim in dims:
        for eps in epss:
            for q in qs:
                dist_fname = "mcgill_orig/dist_matrix_mcgill_original_no_parallel_q_{}_dim_{}.txt".format(q, dim)
                dist_timing_fname = "mcgill_orig/dist_timing_matrix_mcgill_original_q_{}_dim_{}_delta_{}.txt".format(q, dim, eps)
                fname_base = os.path.basename(dist_fname)
                log_name="log_mcgill_timing_orig_correct_spanner_{}_eps_{}.txt".format(fname_base, eps)
                if not os.path.isfile(log_name):
                    continue
                if not os.path.isfile(dist_timing_fname):
                    continue
                bruteforce_time = sum_file(dist_timing_fname)
                print(dist_fname, dist_timing_fname, log_name)
                with open(log_name) as f:
                    for s in f:
                        sr = re.search("blind_greedy_time = ([0-9]+) microseconds", s)
                        if sr:
                            blind_greedy_time = floor( int(sr.group(1)) / 1e6 )
                            final_table[(dim, eps, q)] = (bruteforce_time, blind_greedy_time)


    for dim in dims:
        for eps in epss:
            s = "$\eps = {}$    & ".format(eps)
            for q in qs:
                try:
                    bt, bgt = final_table[(dim, eps, q)]
                except KeyError:
                    bt = bgt = 0
                if q < 3:
                    s += "  {}  &  {}  &&".format(bt, bgt)
                elif q == 3:
                    s += "  {}  &  {}\\\\".format(bt, bgt)
            print(s)
