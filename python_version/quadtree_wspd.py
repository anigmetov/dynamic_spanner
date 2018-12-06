#!/usr/bin/env python3

import numpy as np
import utils_euclidean as ue

dist_matrix = np.zeros(shape=(0, 0))
dist_calls_matrix = np.zeros(shape=(0, 0), dtype=np.bool)

# euclidean_dist_calls = 0
# euclidean_dist_dict = {}

class Point:
    def __init__(self, x, y):
        self.x = x
        self.y = y

def euclidean_dist_with_cache(x, y):
    assert(isinstance(x, int))
    assert(isinstance(y, int))
    # global euclidean_dist_calls
    # global euclidean_dist_dict
    global dist_matrix
    global dist_calls_matrix
    if x == y:
        return 0.0
    dist_calls_matrix[x][y] = dist_calls_matrix[y][x] = 1
    return dist_matrix[x][y]

    # key1 = (tuple(x), tuple(y))
    # if key1 in euclidean_dist_dict:
    #     return euclidean_dist_dict[key1]
    # else:
    #     key2 = (tuple(y), tuple(x))
    #     result = np.sqrt(np.sum((x - y) * (x - y)))
    #     euclidean_dist_calls += 1
    #     euclidean_dist_dict[key1] = result
    #     euclidean_dist_dict[key2] = result
    #     return result


def euclidean_dist_no_count(x, y):
    return np.sqrt(np.sum((x - y) * (x - y)))


class QuadtreeNode:

    def __init__(self, min_x, min_y, max_x, max_y, level, points):
        global global_points
        assert (len(points) >= 1)
        self.min_x = min_x
        self.min_y = min_y
        self.max_x = max_x
        self.max_y = max_y
        self.points = points
        self.level = level

        mid_x = 0.5 * (min_x + max_x)
        mid_y = 0.5 * (min_y + max_y)

        if len(points) == 1:
            return

        p00 = []
        p01 = []
        p10 = []
        p11 = []

        for p in points:
            p_x = global_points[p][0]
            p_y = global_points[p][1]
            if p_x < mid_x:
                if p_y < mid_y:
                    p00.append(p)
                else:
                    p01.append(p)
            else:
                if p_y < mid_y:
                    p10.append(p)
                else:
                    p11.append(p)

        self.children = []
   
        if p00:
            child00 = QuadtreeNode(min_x, min_y, mid_x, mid_y, level + 1, p00)
            self.children.append(child00)

        if p01:
            child01 = QuadtreeNode(min_x, mid_y, mid_x, max_y, level + 1, p01)
            self.children.append(child01)

        if p10:
            child10 = QuadtreeNode(mid_x, min_y, max_x, mid_y, level + 1, p10)
            self.children.append(child10)

        if p11:
            child11 = QuadtreeNode(mid_x, mid_y, max_x, max_y, level + 1, p11)
            self.children.append(child11)

    def is_leaf(self):
        return len(self.points) == 1

    def delta(self):
        if len(self.points) <= 1:
            result = 0
        else:
            result = np.linalg.norm([self.max_x - self.min_x, self.max_y - self.min_y])
        return result

    def corner_points(self):
        p00 = np.array([self.min_x, self.min_y])
        p01 = np.array([self.min_x, self.max_y])
        p10 = np.array([self.max_x, self.min_y])
        p11 = np.array([self.max_x, self.max_y])
        return [p00, p01, p10, p11]

    def distance_to_other(self, other_node):
        return min([euclidean_dist_with_cache(p, p_other) for p in self.points for p_other in other_node.points])

    def distance_to_other_econ(self, other_node):
        # to register one call in cache
        euclidean_dist_with_cache(self.children[0], other_node.children[0])
        return min([np.linalg.norm(p - q) for p in self.corner_points() for q in other_node.corner_points])


    def size(self):
        return len(self.points)


    def print(self):
        print(f"level = {self.level},  ({self.min_x}, {self.min_y}) - ({self.max_x}, {self.max_y}), {len(self.points)}")
        for c in self.children:
            c.print()


class Quadtree:
    def __init__(self, points):
        self.root = QuadtreeNode(0.0, 0.0, 1.0, 1.0, 1, points)



class WspdFromQuadtree:
    def __init__(self, eps, quadtree):
        # global global_points
        self.eps = eps
        self.quadtree = quadtree
        self.wspd = self.gen_wspd(quadtree.root, quadtree.root, 0)
        self.n_points = len(quadtree.root.points)

    def gen_wspd(self, u, v, rec_depth):
        res = set()
        if u.is_leaf() and u == v:
            return res
        if u.delta() < v.delta() or u.delta() == v.delta() and id(u) < id(v):
            return self.gen_wspd(v, u, rec_depth)
        eps = self.eps
        if (1 + eps) * u.delta() <= eps * u.distance_to_other(v):
            res.add((u, v))
            return res
        else:
            return set.union(*[self.gen_wspd(u_c, v, rec_depth + 1)
                               for u_c in u.children])

    def _diameter(self, node):
        global global_points
        return max([euclidean_dist_no_count(global_points[p], global_points[q]) for p in node.points for q in node.points])


    def _dist_between_sets(self, node_1, node_2):
        global global_points
        return min([euclidean_dist_no_count(global_points[p],global_points[q]) for p in node_1.points for q in node_2.points])

    def print_size(self):
        for wspd_mem in self.wspd:
            node_1 = wspd_mem[0]
            node_2 = wspd_mem[1]
            print("sizes {} {}".format(node_1.size(), node_2.size()))

    def check(self):
        # all pairs must be covered
        pair_set = set()
        for wspd_mem in self.wspd:
            node_1 = wspd_mem[0]
            node_2 = wspd_mem[1]
            for ch1 in node_1.points[:]:
                for ch2 in node_2.points[:]:
                    chh1, chh2 = min(ch1, ch2), max(ch1, ch2)
                    if chh1 == 2 and chh2 == 9:
                        print(node_1.points, node_2.points)
                    assert(chh1 != chh2)
                    if (chh1, chh2) in pair_set:
                        print(node_1.points, node_2.points, ch1, ch2, chh1, chh2)
                    assert(not (chh1, chh2) in pair_set)
                    pair_set.add((chh1, chh2))

        all_pair_set = set([(x, y) for x in range(self.n_points)
                            for y in range(x+1, self.n_points)])
        assert (len(all_pair_set.difference(pair_set)) == 0)
        print("WSPD covers all pairs - OK")
        # epsilon-separation
        for wspd_mem in self.wspd:
            node_1 = wspd_mem[0]
            node_2 = wspd_mem[1]
            assert (max(self._diameter(node_1),
                        self._diameter(node_2)) <= self.eps *
                    self._dist_between_sets(node_1, node_2))
        print("WSPD is well-separated - OK")


def make_spanner(eps):
    # global euclidean_dist_calls
    # global euclidean_dist_dict
    global global_points
    global dist_matrix
    global dist_calls_matrix
    qt = Quadtree(range(len(global_points)))
    qt.root.print()
    wspd = WspdFromQuadtree(eps, qt)
    return wspd


if __name__ == "__main__":
    n_points = 100
    n_pairs = int(n_points * (n_points - 1 ) / 2)
    np.random.seed(1)
    global_points = ue.get_uniform_points(n_points, dim=2, max_coord=1.0)
    dist_matrix = ue.distance_matrix(global_points)
    dist_calls_matrix = np.zeros(shape=(n_points, n_points), dtype=np.bool)
    # qt = Quadtree(range(len(global_points)))
    eps = 0.1
    wspd = make_spanner(eps)
    # wspd.print_size()
    # wspd.check()
    print(len(wspd.wspd), n_pairs, len(wspd.wspd) / n_pairs)

