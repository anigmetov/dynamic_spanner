#pragma once

#include <utility>
#include <unordered_set>

#include <boost/config.hpp>
#include <boost/functional/hash.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/incremental_components.hpp>

#include "spdlog/spdlog.h"

#include "cover_tree/cover_tree.h"
#include "cover_tree/wspd_node.h"
#include "cover_tree/boost_graph_import.h"

namespace spd = spdlog;

namespace wasser_spanner {

    struct WSPD
    {
        // types
        using Node = WspdNode;
        using NodePair = std::pair<const Node, const Node>;

        struct NodePairHasher
        {
            std::size_t operator()(const NodePair &p) const
            {
                std::size_t result = 0;
                boost::hash_combine(result, hash_value(p.first));
                boost::hash_combine(result, hash_value(p.second));
                return result;
            }
        };
        using NodePairSet = std::unordered_set<NodePair, NodePairHasher>;

        // data

        CoverTree* m_cover_tree;
        int m_num_points;
        double m_epsilon;
        double m_tau;
        DynamicSpannerR* m_spanner;
        NodePairSet m_wspd;

        std::shared_ptr<spd::logger> console{spd::get("console")};

        std::unique_ptr<Graph> m_graph { nullptr };
        std::vector<size_t> m_rank;
        std::vector<size_t> m_parent;
        DisjointSets m_conn_comps;
        std::unordered_map<Edge, double, EdgeHashStruct> m_distance_cache;

        // methods

        WSPD(CoverTree &ct, double eps);
        void gen_wspd(NodePairSet& result, const Node &a, const Node &b, int rec_depth);
        bool is_valid() const;
        size_t size() const;
        void make_spanner();
        double get_spanner_distance(VertexDescriptor i, VertexDescriptor j);
        double get_max_relative_error();
    };
}
