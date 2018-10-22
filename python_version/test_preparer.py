import os.path
import dionysus as dion
import joblib as jl
import numpy as np


class MyDiagram:
    def __init__(self, dgm_points, id, label=None):
        self.dgm_points = dgm_points
        self.id = id
        self.label = label


def read_diagram(fname, dgm_id):
    pts = []
    with open(fname, 'r') as dgm_file:
        for line in dgm_file:
            x, y = [float(x) for x in line.split()]
            pts.append((x, y))
    return MyDiagram(pts[:], dgm_id, fname)


def get_distance(points_a, points_b, q, delta):
    dgm_a = dion.Diagram(points_a[:])
    dgm_b = dion.Diagram(points_b[:])
    q = int(q)
    if q < 1000:
        return dion.wasserstein_distance(dgm_a, dgm_b, q, delta)
    else:
        return dion.bottleneck_distance(dgm_a, dgm_b, delta)


def read_diagrams_from_dir(dir_name):
    diag_files = [os.path.join(dir_name, f) for f in os.listdir(dir_name) if f.endswith('txt')]
    result = []
    for i in range(len(diag_files)):
        result.append(read_diagram(diag_files[i], i))
    return result


def dist_helper(points_a, points_b, q, delta, i, j):
    return (get_distance(points_a, points_b, q, delta), i, j)


def prepare_distance_matrix(dgms, q, delta):
    n_dgms = len(dgms)
    result = np.zeros((n_dgms, n_dgms), dtype=np.float64)
    res_list = jl.Parallel(jl.delayed(dist_helper)(dgms[i].dgm_points, dgms[j].dgm_points, q, delta, i, j)
                           for i in range(n_dgms)
                           for j in range(i + 1, n_dgms))
    for d, i, j in res_list:
        result[i][j] = result[j][i] = d
    return result
