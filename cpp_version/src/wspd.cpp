#include "cover_tree/wspd.h"

namespace wasser_spanner {

    WSPD::WSPD(CoverTree& ct, double eps) :
            m_cover_tree(&ct),
            m_num_points(ct.m_num_points),
            m_epsilon(eps),
            m_tau(ct.m_base),
            m_spanner(&ct.m_dspanner),
            m_rank(ct.m_num_points),
            m_parent(ct.m_num_points),
            m_conn_comps(&m_rank[0], &m_parent[0])
    {
        Node root_init { m_cover_tree->get_root() };
        gen_wspd(m_wspd, root_init, root_init, 0);
    }

    void WSPD::gen_wspd(NodePairSet& result, const Node& a, const Node& b, int rec_depth)
    {
//        std::string prefix = "[%H:%M:%S.%e]";
//        for(int i = 0; i < rec_depth; ++i) prefix += "  ";
//        prefix += "%v";
//        spd::set_pattern(prefix);

        if (console) { console->debug("Entered gen_wspd, a = {}, b = {}, rec_depth = {}", a, b, rec_depth); }

        if (a.is_leaf_node() and a.get_node_ptr() == b.get_node_ptr()) {
            if (console) { console->info("a is leaf, a == b, return empty set. Exit"); }
            return;
        }

        double a_delta = a.delta();
        double b_delta = b.delta();
        bool swap_needed = a_delta < b_delta or (a_delta == b_delta and a.get_point_idx() > b.get_point_idx());
        if (swap_needed) {
            if (console) {
                console->debug("swap needed, a_delta = {}, b_delta = {}. Recursing, same depth", a_delta, b_delta);
            }
            gen_wspd(result, b, a, rec_depth);
            return;
        }

        assert(a_delta >= b_delta);

        double threshold = (1.0 + m_epsilon) * a_delta / m_epsilon;
        if (m_spanner->is_distance_greater(a.get_point_idx(), b.get_point_idx(), threshold, true)) {
            // TODO - remove debug code
//            auto a_desc = a.get_all_descendants();
//            auto b_desc = b.get_all_descendants();
//            if (a_desc.count(2) == 1 and b_desc.count(3) == 1)
//                console->debug("Rec_depth = {}, adding pair with 2, 3, a: {}, b: {}", rec_depth, a, b);
//
//            if (a_desc.count(3) == 1 and b_desc.count(2) == 1)
//                console->debug("Rec_depth = {}, adding pair with 3, 2, a: {}, b: {}", rec_depth, a, b);
            if (console) { console->debug("Nodes are well-separated, add to result."); }
            result.emplace(a, b);
            return;
        } else {
            if (console) { console->debug("Nodes are not well-separated, recurse to a's children."); }
            for (const auto& a_child : a.get_next_layer_children()) {
                // avoid unnecessary calls
                if (a_child.is_leaf_node() and a_child.get_node_ptr() == b.get_node_ptr()) {
                    continue;
                }
                double a_child_delta = a_child.delta();
                bool swap_needed = a_child_delta < b_delta or
                                   (a_child_delta == b_delta and a_child.get_point_idx() > b.get_point_idx());
                if (swap_needed) {
                    gen_wspd(result, b, a_child, rec_depth + 1);
                } else {
                    gen_wspd(result, a_child, b, rec_depth + 1);
                }
            }
            return;
        }
    }

    bool WSPD::is_valid() const
    {
        // pairs must be well-separated
        for (const auto& node_pair : m_wspd) {
            const Node& node_a = node_pair.first;
            const Node& node_b = node_pair.second;
            if (std::max(node_a.diameter(), node_b.diameter()) > m_epsilon * node_a.get_distance_to_node(node_b)) {
                if (console) { console->critical("not well separated!"); }
                return false;
            }
        }
        if (console) { console->debug("WSPD check: well-separated OK"); }

        // all  pairs must be covered
        std::vector<std::pair<std::unordered_set<int>, std::unordered_set<int>>> unrolled_wspd;

        for (const NodePair& p : m_wspd) {
            unrolled_wspd.emplace_back(p.first.get_all_descendants(), p.second.get_all_descendants());
        }

//        std::cout << "\n----------------------------------------\n";
//        for (const auto& set_pair : unrolled_wspd) {
//            std::cout << "\n{ {";
//            for (const auto p_a : set_pair.first) {
//                std::cout << p_a << ", ";
//            }
//            std::cout << "},  {";
//            for (const auto p_b : set_pair.second) {
//                std::cout << p_b << ", ";
//            }
//            std::cout << "} }";
//        }
//        std::cout << "\n----------------------------------------\n";

        for (int i = 0; i < m_num_points; ++i) {
            for (int j = i + 1; j < m_num_points; ++j) {
                int count = 0;
                for (const auto& set_pair : unrolled_wspd) {
                    if (set_pair.first.count(i) == 1 and set_pair.second.count(j) == 1) {
                        ++count;
                    }
                    if (set_pair.first.count(j) == 1 and set_pair.second.count(i) == 1) {
                        ++count;
                    }
                }
                if (count != 1) {
                    if (console) {
                        console->critical("WSPD is not a decomposition, i = {}, j = {}, count = {}", i, j, count);
                    }
                    return false;
                }
            }
        }
        return true;
    }

    size_t WSPD::size() const
    {
        return m_wspd.size();
    }

    void WSPD::make_spanner()
    {
        m_graph = std::unique_ptr<Graph>(new Graph(m_num_points));
        initialize_incremental_components(*m_graph, m_conn_comps);
        incremental_components(*m_graph, m_conn_comps);
        for (const auto& wspd_node : m_wspd) {
            VertexDescriptor i = wspd_node.first.get_point_idx();
            VertexDescriptor j = wspd_node.second.get_point_idx();
            auto d = m_spanner->get_distance(i, j);
            boost::add_edge(i, j, EdgeWeightProperty(d), *m_graph);
            m_conn_comps.union_set(i, j);
        }
    }

    double WSPD::get_spanner_distance(VertexDescriptor i, VertexDescriptor j)
    {
        if (i == j) {
            return 0.0;
        }
        auto find_res = m_distance_cache.find(std::make_pair(i, j));
        if (find_res != m_distance_cache.end()) {
            return find_res->second;
        } else {
            VertexDescriptor source = boost::vertex(i, *m_graph);
            std::vector<VertexDescriptor> pred_map(m_num_points);
            std::vector<double> dist_map(m_num_points);
            dijkstra_shortest_paths(*m_graph, source,
                                    boost::predecessor_map(&pred_map[0]).distance_map(&dist_map[0]));
            double result = dist_map[j];
            Edge e1 = std::make_pair(i, j);
            Edge e2 = std::make_pair(j, i);
            m_distance_cache[e1] = result;
            m_distance_cache[e2] = result;
            return result;
        }
    }

    double WSPD::get_max_relative_error()
    {
        if (m_graph == nullptr)
            make_spanner();
        double max_relative_error = 0.0;
        for(int i = 0; i < m_num_points; ++i) {
            for (int j = i + 1; j < m_num_points; ++j) {
                double approx_dist = get_spanner_distance(i, j);
                double exact_dist = m_spanner->get_distance(i, j);
                max_relative_error = std::max(max_relative_error, fabs(approx_dist - exact_dist) / exact_dist);
            }
        }
        return max_relative_error;
    }

}
