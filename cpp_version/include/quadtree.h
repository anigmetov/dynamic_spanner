//
// Created by narn on 03.12.18.
//

#ifndef SPANNER_WASSERSTEIN_QUADTREE_H
#define SPANNER_WASSERSTEIN_QUADTREE_H

#include <cmath>
#include <vector>
#include <memory>

namespace quadtree {

using Matrix = std::vector<std::vector<double>>;

struct Point2D {
    double x;
    double y;
};


double dist_2d(const Point2D& p, const Point2D& q)
{
    return std::sqrt( (p.x - q.x) * (p.x - q.x) + (p.y - q.y) * (p.y - q.y) );
}

struct QuadtreeNode {

    using Point = Point2D;

    Point m_min_point;
    Point m_max_point;

    std::vector<Point>& m_all_points;
    std::vector<int>& m_node_point_indices;

    std::unique_ptr<QuadtreeNode> m_child_1 { nullptr };
    std::unique_ptr<QuadtreeNode> m_child_2 { nullptr };
    std::unique_ptr<QuadtreeNode> m_child_3 { nullptr };
    std::unique_ptr<QuadtreeNode> m_child_4 { nullptr };

    std::vector<int> m_children;
    int m_level;

    QuadtreeNode(std::vector<Point>& _points,
                 const std::vector<int>& _point_indices);
};


QuadtreeNode::QuadtreeNode(std::vector<quadtree::QuadtreeNode::Point>& _points,
        const std::vector<int>& _point_indices,
        Point2D p_min,
        Point2D p_max) :
        m_min_point(p_min),
        m_max_point(p_max),

{

}

struct Quadtree {
  Quadtree(const Matrix& distance_matrix, Matrix& distance_computed_matrix);
};

}

#endif //SPANNER_WASSERSTEIN_QUADTREE_H
