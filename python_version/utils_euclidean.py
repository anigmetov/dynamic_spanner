import numpy as np

def euclidean_distance(a, b):
    return np.linalg.norm(a - b)

def get_points(n_points, dim):
    return np.random.randn(n_points, dim)


def get_uniform_points(n_points, dim, max_coord=1.0):
    return np.random.uniform(0.0, max_coord, (n_points, dim))


def get_exponential_points(n_points, dim):
    max_power = 25  # max(n_points / 2.0, 10.0)
    # print(f"max_power = {max_power}")
    shape = (n_points, dim)
    res = np.random.uniform(1.0, max_power, shape)
    return np.power(2.0, res)


def distance_matrix(points):
    n_points = len(points)
    result = np.zeros((n_points, n_points))
    for i in range(n_points):
        for j in range(n_points):
            if i <= j:
                continue
            result[i][j] = result[j][i] = euclidean_distance(points[i, :], points[j, :])
    return result


def distances_list(points):
    n_points = len(points)
    dm = distance_matrix(points)
    return [(dm[i][j], i, j) for i in range(n_points) for j in range(i + 1, n_points)]


def sorted_distances(points):
    return sorted(distances_list(points))


