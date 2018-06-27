#include "cover_tree/wspd_node.h"

namespace wasser_spanner {

// WspdNode

    DynamicSpannerR* WspdNode::dspanner = nullptr;

    WspdNode::WspdNode(CoverTree::Node* node) :
            m_ct_node(node),
            m_level(node->get_level()),
            m_point_idx(node->get_point_idx())
    {
    }

    WspdNode::WspdNode(CoverTree::Node* node, int level) :
            m_ct_node(node),
            m_level(level),
            m_point_idx(node->get_point_idx())
    {
        assert(level <= node->get_level());
    }


    void WspdNode::get_all_descendants_indices_helper(std::unordered_set<int>& res) const
    {
        res.insert(m_point_idx);
        for (const auto& lev_ch_pair : m_ct_node->m_children) {
            if (lev_ch_pair.first < get_level()) {
                for (const auto& ch : lev_ch_pair.second) {
                    WspdNode(ch.get()).get_all_descendants_indices_helper(res);
                }
            }
        }
    }

    std::unordered_set<int> WspdNode::get_all_descendants() const
    {
        std::unordered_set<int> res;
        get_all_descendants_indices_helper(res);
        return res;
    }

    double WspdNode::diameter() const
    {
        double res = 0.0;
        auto indices = get_all_descendants();
        for (auto p_1 = indices.cbegin(); p_1 != indices.cend(); ++p_1) {
            for (auto p_2 = std::next(p_1); p_2 != indices.cend(); ++p_2) {
                res = std::max(res, dspanner->get_distance_no_cache(*p_1, *p_2));
            }
        }
        return res;
    }

    std::vector<WspdNode> WspdNode::get_next_layer_children() const
    {
        std::vector<WspdNode> result;
        result.push_back(get_self_as_child());
        auto find_res = m_ct_node->m_children.find(get_level() - 1);
        if (find_res != m_ct_node->m_children.cend()) {
            for (const auto& c : find_res->second) {
                result.push_back(c.get());
            }
        }
        return result;
    }

    WspdNode WspdNode::get_self_as_child() const
    {
        return WspdNode(m_ct_node, get_level() - 1);
    }

    bool WspdNode::is_leaf_node() const
    {
        for (const auto& lev_ch_pair : m_ct_node->m_children) {
            if (lev_ch_pair.first < get_level()) {
                return false;
            }
        }
        return true;
    }

    double WspdNode::get_distance_to_node(const WspdNode& other) const
    {
        double res = std::numeric_limits<double>::max();
        auto node_a_points = get_all_descendants();
        auto node_b_points = other.get_all_descendants();
        for (const auto& p_a : node_a_points) {
            for (const auto& p_b : node_b_points) {
                res = std::min(res, dspanner->get_distance_no_cache(p_a, p_b));
            }
        }
        return res;
    }

    double WspdNode::delta() const
    {
        if (is_leaf_node()) {
            return 0.0;
        } else {
            return std::pow(base, get_level() + 2);
        }
    }

    std::size_t hash_value(const WspdNode& node)
    {
        std::size_t result = 0;
        boost::hash<int> hasher;
        boost::hash<CoverTree::Node*> hasher_;
        boost::hash_combine(result, hasher(node.get_level()));
        boost::hash_combine(result, hasher(node.get_point_idx()));
        boost::hash_combine(result, hasher_(node.get_node_ptr()));
        return result;
    }

    bool operator==(const WspdNode& a, const WspdNode& b)
    {
        // if two nodes point to the same cover tree node, they have the same point
        assert((a.get_node_ptr() != b.get_node_ptr()) xor (a.get_point_idx() == b.get_point_idx()));
        return a.get_node_ptr() == b.get_node_ptr() and
               a.get_level() == b.get_level();
    }

    std::ostream& operator<<(std::ostream& s, const WspdNode& n)
    {
        s << "WspdNode(point = " << n.get_point_idx() << ", level = " << n.get_level() << ", original level = "
          << n.get_node_ptr()->get_level() << ", is_leaf = " << n.is_leaf_node() << ")";
        return s;
    }
}
