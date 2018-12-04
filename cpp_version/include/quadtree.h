//
// Created by narn on 03.12.18.
//

#ifndef SPANNER_WASSERSTEIN_QUADTREE_H
#define SPANNER_WASSERSTEIN_QUADTREE_H

#include <vector>

namespace quadtree {

using Matrix = std::vector<std::vector<double>>;

struct QuadtreeNode {
  std::vector<int> m_children;
  int m_level;
};


struct Quadtree {
  Quadtree(const Matrix& distance_matrix, Matrix& distance_computed_matrix);
};



}

#endif //SPANNER_WASSERSTEIN_QUADTREE_H
