#ifndef SPANNER_WASSERSTEIN_WSPD_NODE_H
#define SPANNER_WASSERSTEIN_WSPD_NODE_H

#include <iostream>
#include <vector>
#include <memory>
#include <unordered_set>
#include <unordered_map>

#include <boost/functional/hash.hpp>

#include "cover_tree/cover_tree.h"

namespace wasser_spanner {

    class WspdNode
    {
    private:
        // members
        CoverTree::Node* m_ct_node { nullptr };
        const int m_level { CoverTree::k_max_lvl };
        const int m_point_idx { CoverTree::k_invalid_idx };

        std::unordered_map<int, CoverTree::NodeVec>* m_children { nullptr };

        void get_all_descendants_indices_helper(std::unordered_set<int>& res) const;

        WspdNode get_self_as_child() const;

    public:
        static DynamicSpannerR* dspanner;
        static double constexpr base { 2.0 };

        WspdNode(CoverTree::Node* node);

        WspdNode(CoverTree::Node* node, int level);

        int get_level() const
        { return m_level; }

        int get_point_idx() const
        { return m_point_idx; }

        CoverTree::Node* get_node_ptr() const
        { return m_ct_node; }

        bool is_leaf_node() const;

        std::unordered_set<int> get_all_descendants() const;

        double delta() const;

        double get_distance_to_node(const WspdNode& other) const;

        double diameter() const;

        std::vector<WspdNode> get_next_layer_children() const;
    };

    std::size_t hash_value(const WspdNode& node);

    bool operator==(const WspdNode& a, const WspdNode& b);

    std::ostream& operator<<(std::ostream& s, const WspdNode& n);

} // namespace wasser_spanner
#endif //SPANNER_WASSERSTEIN_WSPD_NODE_H
