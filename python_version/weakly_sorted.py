#/usr/bin/env python3

import numpy as np
import random

class QuasiSortedDistances:
    def __init__(self, qsorted_distances, n_buckets):
        self.qsorted_distances = qsorted_distances
        self.n_buckets = n_buckets


def get_buckets(min_dist, max_dist):
    n_buckets = int(np.ceil(np.log2(max_dist / min_dist)))
    scaling_factor = min_dist
    result = [(scaling_factor * (2 ** k), 4.0 * scaling_factor * (2 ** k)) for k in range(n_buckets)]
    return result


def get_bucket_index(buckets, d):
    dist = d[0]
    cand_list = []
    for i in range(len(buckets)):
        lb, ub = buckets[i]
        if dist >= lb and dist <= ub:
            cand_list.append(i)
    if len(cand_list) == 0:
        b0 = buckets[0]
        b1 = buckets[-1]
        print(f"ACHTUNG! dist = {dist}, {b0} - {b1}  bucket not found")
    assert (len(cand_list) > 0)
    if len(cand_list) == 1:
        return cand_list[0]
    else:
        idx = random.randint(0, len(cand_list) - 1)
        return cand_list[idx]


def quasi_sorted_distances(points, dist_list=None):
    if dist_list is None:
        n_points = len(points)
        dm = distance_matrix(points)
        dist_list = [(dm[i][j], i, j) for i in range(n_points) for j in range(i + 1, n_points)]
    max_dist = max(dist_list)[0]
    min_dist = min(dist_list)[0]
    buckets = get_buckets(min_dist, max_dist)
    n_buckets = len(buckets)
    result = [[] for i in range(n_buckets)]
    for d in dist_list:
        result[get_bucket_index(buckets, d)].append(d)
    map(random.shuffle, result)
    return QuasiSortedDistances([x for y in result for x in y], n_buckets)


