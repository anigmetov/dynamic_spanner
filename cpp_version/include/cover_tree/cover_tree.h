#pragma once

#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <unordered_set>

#include "spdlog/spdlog.h"
#include "wasserstein/wasserstein_space_point.h"

namespace spd = spdlog;

namespace wasser_spanner {

    using DynamicSpannerR = DynamicSpanner<double>;
    using AuctionParamsR = hera::AuctionParams<double>;
    using Point = DynamicSpannerR::WassersteinSpacePoint;
    using MatrixR = DynamicSpannerR::MatrixR;


    struct CoverTree
    {
        static constexpr int k_max_lvl = std::numeric_limits<int>::max() - 10;
        static constexpr int k_invalid_idx = -1;

        struct Node;

        using NodeVec = std::vector<Node*>;

        struct Node
        {
            // data
            int m_level { k_max_lvl + 1 };
            int m_point_idx { k_invalid_idx };
            // m_children vector does not include the node itself, it is implicit
            std::unordered_map<int, std::vector<std::unique_ptr<Node>>> m_children;

            // methods
            Node(int _point_idx, int _level);

            int get_point_idx() const
            { return m_point_idx; }

            int get_level() const
            { return m_level; }

            int get_num_children() const;

            NodeVec get_children(int level) const;

            friend std::ostream& operator<<(std::ostream& os, const Node& node);

            std::shared_ptr<spd::logger> console { spd::get("console") };
        };

        // data members
        static constexpr double m_base { 2.0 };
        std::vector<Point> m_points;
        unsigned int m_num_points;
        int m_max_level;
        int m_min_level;
        std::unique_ptr<Node> m_root;
        DynamicSpannerR m_dspanner;
        unsigned int m_num_nodes { 0 };
        std::unordered_map<int, Node*> m_point_to_node;

        std::shared_ptr<spd::logger> console { spd::get("console") };


        // member functions
        CoverTree(double max_dist,
                  const std::vector<Point>& points,
                  const MatrixR& distance_matrix = MatrixR(),
                  const AuctionParamsR& _ap = AuctionParamsR());

        void insert(const Point& new_point);

        bool insert_rec(const Point& p,
                        const std::vector<Node*>& Qi,
                        const int level,
                        bool& inserted);

        void add_child(Node* node, int level, int point_idx);

        Node* get_root() const;

        // auxiliary functions
        bool are_all_children_further(const Point& p,
                                      const NodeVec& Qi,
                                      double threshold,
                                      int level);

        NodeVec get_near_points(const Point& p, const NodeVec& Q, double threshold, int level);

        bool is_some_node_closer(const Point& p,
                                 const NodeVec& near_points,
                                 double threshold,
                                 Node*& cand_parent);

        bool is_valid_tree() const;

        bool is_valid_node(const Node* node, std::set<int>& remaining_points) const;

        void separation_invariant_helper(const Node* node, std::unordered_map<int, NodeVec>& all_nodes) const;

        bool is_separation_invariant_satisfied() const;

        void print_node(std::ostream& s, CoverTree::Node* node, int indent) const;

        void print_tree(std::ostream& s) const;
    };

    // CoverTree class


    bool operator==(const CoverTree::Node& a, const CoverTree::Node& b);
}
