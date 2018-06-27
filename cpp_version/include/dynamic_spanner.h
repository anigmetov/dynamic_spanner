#ifndef DYNAMIC_SPANNER_H
#define DYNAMIC_SPANNER_H
/*
#include <vector>
#include <utility>

#include <boost/config.hpp>
#include <boost/functional/hash.hpp>

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/incremental_components.hpp>

using namespace boost;

// class Point must provide function distance
template<class Point, class Real = double>
class DynamicSpanner
{
public:

    using EdgeWeightProperty = property<edge_weight_t, double>;
    using Graph = adjacency_list<vecS, vecS,undirectedS, no_property, EdgeWeightProperty>;
    using VertexDescriptor = typename graph_traits<Graph>::vertex_descriptor;
    using EdgeDescriptor = typename graph_traits<Graph>::edge_descriptor;
    using EdgeIterator = typename graph_traits<Graph>::edge_iterator;
    using Edge = typename std::pair<VertexDescriptor, VertexDescriptor>;
    using DisjointSets = disjoint_sets<size_t*, size_t*>;
    using WeightMap = property_map<Graph, edge_weight_t>::type;

    struct EdgeHashStruct {
        size_t operator()(const Edge& p) const
        {
            size_t seed = 0;
            boost::hash_combine(seed, p.first);
            boost::hash_combine(seed, p.second);
            return seed;
        }
    };

    DynamicSpanner(std::vector<Point>& _pts) :
        points_(_pts),
        num_points_(_pts.size()),
        graph_(_pts.size()),
        rank_(num_points_),
        parent_(num_points_),
        conn_comps_(&rank_[0], &parent_[0])
    {
        initialize_incremental_components(graph_, conn_comps_);
        incremental_components(graph_, conn_comps_);
    }

    std::pair<Real, Real> distance_estimate(VertexDescriptor i, VertexDescriptor j)
    {
        Real d;
        if (is_distance_cached(i, j, d)) {
            return std::make_pair(d, d);
        } else {
            return std::make_pair(lower_bound(i,j), upper_bound(i,j));
        }
    }

    // upper bound for distance between points_[i] and points_[j]
    Real upper_bound(VertexDescriptor i, VertexDescriptor j)
    {
        VertexDescriptor source = vertex(i, graph_);
        std::vector<VertexDescriptor> pred_map(num_points_);
        std::vector<Real> dist_map(num_points_);
        dijkstra_shortest_paths(graph_, source,
                predecessor_map(&pred_map[0]).distance_map(&dist_map[0]));
        Real result = dist_map[j];
        return result;
    }

    bool is_distance_less(VertexDescriptor i, VertexDescriptor j,
            Real value, bool strict)
    {
        Real d;
        if (is_distance_cached(i, j, d)) {
            return d < value or (not strict and d <= value);
        }

        Real ub = upper_bound(i, j);
        if (ub < value or (not strict and ub <= value))
            return true;
        Real lb = lower_bound(i, j);
        if (lb >= value or (strict and lb > value))
            return false;
        d = points_[i].distance(points_[j]);
        add_edge(i, j, d);
        return d < value or (not strict and d <= value);
    }

    bool is_distance_greater(VertexDescriptor i, VertexDescriptor j,
            Real value, bool strict)
    {
        return not is_distance_less(i, j, value, not strict);
    }


    bool is_distance_cached(VertexDescriptor i, VertexDescriptor j, Real& d) const
    {
        auto find_res = distance_cache_.find(std::make_pair(i,j));
        if (find_res != distance_cache_.end()) {
            d = find_res->second;
            return true;
        } else {
            return false;
        }
    }

    bool is_distance_cached(VertexDescriptor i, VertexDescriptor j) const
    {
        Real d;
        return is_distance_cached(i, j, d);
    }


    Real lower_bound(VertexDescriptor i, VertexDescriptor j)
    {
        double result = 0.0;

        EdgeIterator ei, ei_end;
        WeightMap weight_map = get(edge_weight, graph_);

        for(boost::tie(ei, ei_end) = edges(graph_); ei != ei_end; ++ei) {
            VertexDescriptor u = source(*ei, graph_);
            VertexDescriptor v = target(*ei, graph_);
            if (same_component(u, i, conn_comps_)) {
                double ui_dist = upper_bound(u, i, graph_);
                double vj_dist = upper_bound(v, j, graph_);
                double uj_dist = upper_bound(u, j, graph_);
                double vi_dist = upper_bound(v, i, graph_);
                double uv_dist = get(weight_map, *ei);
                result = std::max(result, uv_dist - ui_dist - vj_dist);
                result = std::max(result, uv_dist - uj_dist - vi_dist);
            }
        }
        return result;
    }

    void add_edge(VertexDescriptor i, VertexDescriptor j, Real distance)
    {
        add_edge(i, j, EdgeWeightProperty(distance), graph_);
        conn_comps_.union_set(i, j);
        Edge edge_1 = std::make_pair(i, j);
        Edge edge_2 = std::make_pair(j, i);
        distance_cache_[edge_1] = distance;
        distance_cache_[edge_2] = distance;
    }

private:
    std::vector<Point> points_;
    size_t num_points_;
    Graph graph_;

    // for maintaining connected components we need rank_ and parent_
    std::vector<size_t> rank_;
    std::vector<size_t> parent_;
    DisjointSets conn_comps_;

    std::unordered_map<Edge, Real, EdgeHashStruct> distance_cache_;

};
*/
#endif
