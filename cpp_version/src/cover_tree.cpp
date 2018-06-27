#include <boost/functional/hash.hpp>

#include "cover_tree/cover_tree.h"

namespace wasser_spanner {

    CoverTree::Node::Node(int _point_idx, int _level) :
            m_level(_level),
            m_point_idx(_point_idx)
    {
        if (console) { console->debug("Created Node, level = {}, point_idx = {}", m_level, m_point_idx); }
    }

    bool operator==(const CoverTree::Node& a, const CoverTree::Node& b)
    {
        assert(a.m_point_idx != b.m_point_idx or (a.m_level == b.m_level and a.m_children == b.m_children));
        return a.m_level == b.m_level and a.m_point_idx == b.m_point_idx;
    }

    int CoverTree::Node::get_num_children() const
    {
        int result = 0;
        for (const auto& lev_ch_pair : m_children)
            result += lev_ch_pair.second.size();
        return result;
    }

    void CoverTree::add_child(Node* node, int level, int point_idx)
    {
        if (console) {
            console->debug("Enter add child, level = {}, point_idx = {}, parent node = {}", level, point_idx,
                           node->m_point_idx);
        }
        assert(m_point_to_node.find(point_idx) == m_point_to_node.end());
        assert(node->get_level() > level);

        node->m_children[level].emplace_back(new Node(point_idx, level));
        m_point_to_node[point_idx] = node->m_children[level].back().get();
    }


    CoverTree::NodeVec CoverTree::Node::get_children(int level) const
    {
        std::vector<CoverTree::Node*> result;

        auto find_res = m_children.find(level);
        if (find_res != m_children.end()) {
            for (const auto& c : find_res->second) {
                result.push_back(c.get());
            }
        }

        return result;
    }

    bool CoverTree::are_all_children_further(const Point& p,
                                             const CoverTree::NodeVec& Qi,
                                             double threshold,
                                             int level)
    {
        if (console) { console->debug("Enter are_all_children_further, id = {}, threshold = {}", p.id, threshold); }

        bool use_lazy = true;

        if (use_lazy) {
            bool lazy_failed = false;

            for (const auto node : Qi) {
                if (console) { console->debug("lazy considering node {}", *node); }

                // node with point p has point p as a child - check separately
                int lazy_comp_result = m_dspanner.is_distance_less_wo_hera(p.id, node->m_point_idx, threshold, false);
                if (1 == lazy_comp_result) {
                    if (console) {
                        console->debug("Exit are_all_children_further, id = {}, threshold = {}, lazy false-1", p.id,
                                       threshold);
                    }
                    return false;
                } else if (-1 == lazy_comp_result) {
                    lazy_failed = true;
                    break;
                }

                if (not lazy_failed) {
                    // check proper children
                    for (const auto& child_node : node->get_children(level)) {

                        int lazy_comp_result = m_dspanner.is_distance_less_wo_hera(p.id, child_node->m_point_idx,
                                                                                   threshold, false);
                        if (1 == lazy_comp_result) {
                            if (console) {
                                console->debug("Exit are_all_children_further, id = {}, threshold = {}, lazy false-2",
                                               p.id, threshold);
                            }
                            return false;
                        } else if (-1 == lazy_comp_result) {
                            lazy_failed = true;
                            break;
                        }
                    }
                }
            }

            if (not lazy_failed) {
                if (console) {
                    console->debug("Exit are_all_children_further, id = {}, threshold = {}, lazy true", p.id,
                                   threshold);
                }
                return true;
            }
        }

        for (const auto node : Qi) {
            if (console) { console->debug("non-lazy considering node {}", *node); }

            if (m_dspanner.is_distance_less(p.id, node->m_point_idx, threshold, false)) {
                if (console) {
                    console->debug("Exit are_all_children_further, id = {}, threshold = {}, false-1", p.id, threshold);
                }
                return false;
            }

            for (const auto child_node : node->get_children(level)) {

                if (console) { console->debug("In are_all_children_further, processing child {}", *child_node); }

                if (m_dspanner.is_distance_less(p.id, child_node->m_point_idx, threshold, false)) {
                    if (console) {
                        console->debug("Exit are_all_children_further, id = {}, threshold = {}, false-2", p.id,
                                       threshold);
                    }
                    return false;
                }
            }
        }

        if (console) {
            console->debug("Exit are_all_children_further, id = {}, threshold = {}, true", p.id, threshold);
        }
        return true;
    }

    CoverTree::NodeVec CoverTree::get_near_points(const Point& p,
                                                  const CoverTree::NodeVec& Q,
                                                  double threshold,
                                                  int level)
    {
        CoverTree::NodeVec result;
        for (const auto& node : Q) {

            // node is its own child - separate check
            if (m_dspanner.is_distance_less(p.id, node->m_point_idx, threshold, false)) {
                result.push_back(node);
            }

            // loop over proper children
            for (const auto& child_node : node->get_children(level)) {
                if (m_dspanner.is_distance_less(p.id, child_node->m_point_idx, threshold, false)) {
                    result.push_back(child_node);
                }
            }
        }

        return result;
    }

    bool CoverTree::is_some_node_closer(const Point& p,
                                        const CoverTree::NodeVec& near_points,
                                        double threshold,
                                        Node*& cand_parent)
    {
        for (const auto& node : near_points) {
            if (m_dspanner.is_distance_less(p.id, node->m_point_idx, threshold, false)) {
                cand_parent = node;
                return true;
            }
        }
        return false;
    }

    bool CoverTree::insert_rec(const Point& p,
                               const std::vector<Node*>& Qi,
                               const int level,
                               bool& inserted)
    {
        if (console) {
            console->debug("Calling recursive insert, id = {}, level = {}, Qi:", p.id, level);
            for (const auto q : Qi) {
                console->debug("Node {}", *q);
            }
        }

        double sep = std::pow(m_base, level);
        if (console) { console->debug("Sep = {}", sep); }

        if (are_all_children_further(p, Qi, sep, level - 1)) {
            if (console) {
                console->debug("Exit recursive insert, id = {}, level = {}, false, all children are further", p.id,
                               level);
            }
            return false;
        }

        if (console) { console->debug("are_all_children_further is false, continue"); }

        auto near_points = get_near_points(p, Qi, sep, level - 1);

//        if (console) {
//            for (const auto np : near_points) {
//                console->debug("Near point is {}, distance {}", *np,
//                               dspanner.get_distance_no_cache(np->get_point_idx(), p.id));
//            }
//        }

        Node* cand_parent = nullptr;
        if (not insert_rec(p, near_points, level - 1, inserted) and not inserted and
            is_some_node_closer(p, Qi, sep, cand_parent)) {
            assert(cand_parent != nullptr);
            if (console) {
                console->debug("cand_parent is {}", *cand_parent);
            }
            add_child(cand_parent, level - 1, p.id);
            inserted = true;
            m_min_level = std::min(m_min_level, level - 1);
            m_num_nodes++;
            if (console) { console->debug("Exit recursive insert, id = {}, level = {}, true", p.id, level); }
            return true;
        } else {
            if (console) { console->debug("Exit recursive insert, id = {}, level = {}, false", p.id, level); }
            return false;
        }
    }

    CoverTree::CoverTree(double max_dist, const std::vector<Point>& _points, const MatrixR& dist_matrix,
                         const AuctionParamsR& _ap) :
            m_points(_points),
            m_num_points(_points.size()),
            m_max_level(std::ceil(std::log(max_dist) / std::log(m_base))),
            m_min_level(m_max_level),
            m_dspanner(m_points, dist_matrix, _ap)
    {
        for (const auto& p : m_points) {
            if (console) { console->debug("in CoverTree ctor, processing point, id = {}", p.id); }
            insert(p);
        }
    }

    void CoverTree::insert(const Point& new_point)
    {
        bool inserted = false;
        if (console) { console->debug("Enter insert, id = {}", new_point.id); }
        if (m_root) {
            if (not insert_rec(new_point, NodeVec(1, m_root.get()), m_max_level, inserted) and not inserted) {
                add_child(m_root.get(), m_max_level - 1, new_point.id);
            }
        } else {
            if (console) { console->debug("Creating root"); }
            m_root = std::unique_ptr<Node>(new Node(new_point.id, m_max_level));
            m_point_to_node[new_point.id] = m_root.get();
            m_num_nodes++;
        }

        if (console) {
            console->debug("Exit insert, id = {}", new_point.id);
            console->debug("Computed distances: {}", m_dspanner.m_hera_distance_calls);
        }
    }

    CoverTree::Node* CoverTree::get_root() const
    {
        return m_root.get();
    }


    void CoverTree::separation_invariant_helper(const CoverTree::Node* node,
                                                std::unordered_map<int, NodeVec>& all_nodes) const
    {
        for (const auto& lev_ch_pair : node->m_children) {
            int level = lev_ch_pair.first;
            assert(level < node->get_level());
            for (const auto& c : lev_ch_pair.second) {
                all_nodes[level].push_back(c.get());
                separation_invariant_helper(c.get(), all_nodes);
            }
        }
    }

    bool CoverTree::is_separation_invariant_satisfied() const
    {
        std::unordered_map<int, NodeVec> all_nodes;
        separation_invariant_helper(get_root(), all_nodes);
        for (const auto& lev_nodes_pair : all_nodes) {
            int level = lev_nodes_pair.first;
            auto& node_vec = lev_nodes_pair.second;
            double sep = std::pow(m_base, level);
            for (size_t i = 0; i < node_vec.size(); ++i) {
                for (size_t j = i + 1; j < node_vec.size(); ++j) {
                    double d = m_dspanner.get_distance_no_cache(node_vec[i]->get_point_idx(),
                                                                node_vec[j]->get_point_idx());
                    if (d <= sep) {
                        if (console) {
                            console->critical(
                                    "Separation invariant failed, level = {}, node = {}, {}, d = {}, sep = {}", level,
                                    node_vec[i]->m_point_idx, node_vec[j]->m_point_idx, d, sep);
                        }
                        return false;
                    }
                }
            }
        }
        return true;
    }

    bool CoverTree::is_valid_node(const CoverTree::Node* node, std::set<int>& remaining_points) const
    {
        if (not(0 <= node->get_point_idx() and node->get_point_idx() < (int) m_num_points)) {
            if (console) { console->critical("node w/o point {}", *node); }
            return false;
        }

        if (node->get_level() > k_max_lvl) {
            if (console) { console->critical("node with invalid level {}", *node); }
            return false;
        }

        const int point_idx = node->m_point_idx;
        remaining_points.erase(point_idx);

        for (const auto& lev_ch_pair : node->m_children) {
            int level = lev_ch_pair.first;
            for (const auto& pchild : lev_ch_pair.second) {
                if (level != pchild->m_level) {
                    if (console) {
                        console->critical(
                                "Level in node and parent's m_children mismatch, level = {}, node = {}, node.m_level = {}",
                                level, *node);
                    }
                }
                assert(level == pchild->m_level);
                int child_idx = pchild->m_point_idx;
                double dist_to_child = m_dspanner.get_distance_no_cache(point_idx, child_idx);
                double sep = std::pow(m_base, level + 1);

                if (dist_to_child > sep) {
                    if (console) { console->critical("covering invariant failed {}", *node); }
                    return false;
                }

                if (not is_valid_node(pchild.get(), remaining_points)) {
                    return false;
                }
            }
        }
        return true;
    }

    bool CoverTree::is_valid_tree() const
    {
        // all m_points must be in the tree
        std::set<int> point_indices;
        for (int i = 0; i < (int) m_num_points; ++i) {
            // hinted insertion
            point_indices.insert(point_indices.end(), i);
        }

        for (int i = 0; i < (int) m_num_points; ++i) {
            auto find_res = m_point_to_node.find(i);
            if (find_res == m_point_to_node.cend()) {
                if (console) { console->critical("cannot find point {}", i); }
                return false;
            };
            if (find_res->second->m_point_idx != i) {
                if (console) { console->critical("error in dict at point {}", i); }
                return false;
            }
        }

        // check covering invariant
        if (not(is_valid_node(get_root(), point_indices) and point_indices.empty())) {
            return false;
        }

        if (console) { console->info("Covering OK"); }

        return is_separation_invariant_satisfied();
    }


    std::ostream& operator<<(std::ostream& os, const CoverTree::Node& node)
    {
        os << " Node(level = " << node.m_level << ", point idx = " << node.m_point_idx << ", #children = "
           << node.get_num_children() << ") ";
        return os;
    }


    void CoverTree::print_tree(std::ostream& s) const
    {
        s << "CoverTree. Num_points = " << m_num_points << ", num_nodes = " << m_num_nodes << " Root: " << std::endl;
        print_node(s, get_root(), 0);
        s << std::endl;
    }

    void CoverTree::print_node(std::ostream& s, CoverTree::Node* node, int indent) const
    {
        s << "\n";
        for (int i = 0; i < indent; ++i)
            s << " ";
        s << *node;
        for (const auto& lev_ch_pair : node->m_children) {
            for (const auto& c : lev_ch_pair.second) {
                print_node(s, c.get(), 4 + indent);
            }
        }
    }

} // namespace wasser_spanner
