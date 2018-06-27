#ifndef WASSERSTEIN_SPACE_POINT_H
#define WASSERSTEIN_SPACE_POINT_H

#include <unordered_map>
#include <iostream>
#include <utility>
#include <memory>

#include <boost/config.hpp>
#include <boost/functional/hash.hpp>

#include "spdlog/spdlog.h"

#include "wasserstein/wasserstein.h"
#include "bottleneck/bottleneck.h"
#include "cover_tree/boost_graph_import.h"

namespace spd = spdlog;

constexpr int INVALID_ID = -1;

// class Point must provide function distance
template<class Real_ = double>
class DynamicSpanner
{
public:

// types
    using Real = Real_;
    using DiagramPointR = std::pair<Real, Real>;
    using DiagramR = std::vector<DiagramPointR>;
    using AuctionParamsR = typename hera::AuctionParams<Real>;
    using MatrixR = typename std::vector<std::vector<Real>>;

    struct WassersteinSpacePoint
    {
        // diagrams are stored without any diagonal m_points
        // Hera will add projections itself
        //
        using Real = Real_;
        using IdType = int;
        //using CacheKeyType = std::pair<IdType, IdType>;
        //using CacheType = std::unordered_map<CacheKeyType, Real, boost::hash<CacheKeyType>>;

        DiagramR dgm;
        IdType id;
        std::string comment;

        DynamicSpanner<Real>* dyn_spanner;

        WassersteinSpacePoint() :
                id(INVALID_ID)
        {
        }


        WassersteinSpacePoint(const DiagramR& _dgm, int _id) :
                dgm(_dgm), id(_id)
        {
        }

        // read diagram from text file
        WassersteinSpacePoint(const std::string& dgm_file, int _id) :
                id(_id),
                comment(dgm_file)
        {
            if (!hera::read_diagram_point_set(dgm_file, dgm)) {
                throw std::runtime_error(std::string("Cannot read file ") + dgm_file);
            }
        }

        std::string to_string() const
        {
            std::string result = "Point(id = ";
            return result;
        }

        void print() const
        {
        };

        bool operator==(const WassersteinSpacePoint& other) const
        {
            return id == other.id;
        };


    };

    using Point = WassersteinSpacePoint;

private:


// data

public:
    bool m_use_geom_lower_bound { true };
    int m_hera_distance_calls { 0 };
    int m_spared_hera_calls { 0 };
    int m_cached_calls { 0 };
    int m_geom_lower_bound_useful { 0 };
    std::vector<Point> m_points;
    size_t m_num_points;
    Graph m_graph;

    // for maintaining connected components we need m_rank and m_parent
    std::vector<size_t> m_rank;
    std::vector<size_t> m_parent;
    DisjointSets m_conn_comps;


    // if we have precomputed distance matrix for experiment,
    // use it; it it's empty, compute on the fly
    MatrixR m_distance_matrix;

    // wasserstein-specific

    AuctionParamsR m_auction_params;

    std::unordered_map<Edge, Real, EdgeHashStruct> m_distance_cache;
    std::unordered_map<Edge, Real, EdgeHashStruct> m_upper_bound_cache;
    std::unordered_map<Edge, Real, EdgeHashStruct> m_lower_bound_cache;
    std::unordered_map<Edge, Real, EdgeHashStruct> m_geom_lower_bound_cache;
    std::unordered_set<Edge, EdgeHashStruct> m_requested_distances;

    std::shared_ptr<spd::logger> console { spd::get("console") };

// methods
public:
    DynamicSpanner(std::vector<Point>& _pts,
                   const MatrixR& distance_matrix,
                   AuctionParamsR _auction_params) :
            m_hera_distance_calls(0),
            m_spared_hera_calls(0),
            m_cached_calls(0),
            m_points(_pts),
            m_num_points(_pts.size()),
            m_graph(_pts.size()),
            m_rank(m_num_points),
            m_parent(m_num_points),
            m_conn_comps(&m_rank[0], &m_parent[0]),
            m_distance_matrix(distance_matrix),
            m_auction_params(_auction_params)
    {
        initialize_incremental_components(m_graph, m_conn_comps);
        incremental_components(m_graph, m_conn_comps);
    }

    // upper bound for distance between m_points[i] and m_points[j]
    Real upper_bound(VertexDescriptor i, VertexDescriptor j)
    {
        Real result;
        if (console) { console->debug("Enter upper_bound, i = {}, j = {}", i, j); }

        if (is_distance_cached(i, j, result)) {
            if (console) { console->debug("Exit upper_bound, i = {}, j = {}, result = {} cached", i, j, result); }
            return result;
        }

        VertexDescriptor source = vertex(i, m_graph);
        std::vector<VertexDescriptor> pred_map(m_num_points);
        std::vector<Real> dist_map(m_num_points);
        dijkstra_shortest_paths(m_graph, source,
                                boost::predecessor_map(&pred_map[0]).distance_map(&dist_map[0]));
        result = dist_map[j];
        if (console) { console->debug("Exit upper_bound, i = {}, j = {}, result = {}", i, j, result); }

        m_upper_bound_cache[std::make_pair(i, j)] = result;
        m_upper_bound_cache[std::make_pair(j, i)] = result;

        return result;
    }

    // call Hera if necessary, do not cache result or change stats
    double get_distance_no_cache(VertexDescriptor i, VertexDescriptor j) const
    {
        if (m_distance_matrix.empty()) {
            DiagramR d1 = m_points[i].dgm;
            DiagramR d2 = m_points[j].dgm;
            if (m_auction_params.wasserstein_power != std::numeric_limits<double>::infinity()
                and m_auction_params.wasserstein_power != hera::get_infinity()) {
                    return hera::wasserstein_dist(d1, d2, m_auction_params);
            } else {
                    return hera::bottleneckDistApprox(d1, d2, m_auction_params.delta);
            }
        } else {
            return m_distance_matrix[i][j];
        }
    }


    double get_distance(VertexDescriptor i, VertexDescriptor j)
    {
        console->debug("get_distance called for {} {}", i, j);
        double result;
        if (is_distance_cached(i, j, result)) {
            m_cached_calls++;
        } else {
            result = get_distance_no_cache(i, j);
            m_hera_distance_calls++;
            cache_distance(i, j, result);
        }
        add_edge(i, j, result);
        return result;
    }

    Real get_cached_upper_bound(VertexDescriptor i, VertexDescriptor j) const
    {
        Edge e = std::make_pair(i, j);
        auto find_res = m_upper_bound_cache.find(e);
        if (find_res != m_upper_bound_cache.cend()) {
            return find_res->second;
        } else {
            return std::numeric_limits<Real>::max();
        }
    }

    Real get_geom_lower_bound(VertexDescriptor i, VertexDescriptor j)
    {
        if (m_auction_params.wasserstein_power == hera::get_infinity()) {
            return 0.0;
        }
        Edge e = std::make_pair(i, j);
        auto find_res = m_geom_lower_bound_cache.find(e);
        if (find_res != m_geom_lower_bound_cache.cend()) {
            return find_res->second;
        } else {
            DiagramR d1 = m_points[i].dgm;
            DiagramR d2 = m_points[j].dgm;
            Real result = hera::wasserstein_dist_lower_bound(d1, d2, m_auction_params);
            if (result > 1.001 * get_distance_no_cache(i, j)) {
                console->critical("geom. lower bound = {}, dist = {}, i = {}, j = {}", result,
                                  get_distance_no_cache(i, j), i, j);
            }
            assert(result <= 1.001 * get_distance_no_cache(i, j));
            Edge e1 = std::make_pair(j, i);
            m_geom_lower_bound_cache[e] = result;
            m_geom_lower_bound_cache[e1] = result;
            return result;
        }
    }

    Real get_cached_lower_bound(VertexDescriptor i, VertexDescriptor j) const
    {
        Edge e = std::make_pair(i, j);
        auto find_res = m_lower_bound_cache.find(e);
        if (find_res != m_lower_bound_cache.cend()) {
            return find_res->second;
        } else {
            return 0.0;
        }
    }

    bool is_distance_less(VertexDescriptor i, VertexDescriptor j,
                          Real value, bool strict)
    {

        if (i == j) {
            if (strict) {
                return 0.0 < value;
            } else {
                return 0.0 <= value;
            }
        }

        m_requested_distances.emplace(i, j);
        m_requested_distances.emplace(j, i);

        if (console) { console->debug("Enter is_distance_less, i = {}, j = {}", i, j); }

        Real d;

        if (is_distance_cached(i, j, d)) {
            m_cached_calls++;
            if (console) { console->debug("Exit is_distance_less, i = {}, j = {}, cached", i, j); }
            return d < value or (not strict and d <= value);
        } else {
            d = 0.0;
        }

        Real ub = get_cached_upper_bound(i, j);
        if (ub < value or (not strict and ub <= value)) {
            m_spared_hera_calls++;
            if (console) { console->debug("Exit is_distance_less, i = {}, j = {}, spared ub", i, j); }
            return true;
        } else {
            ub = upper_bound(i, j);
            if (ub < value or (not strict and ub <= value)) {
                m_spared_hera_calls++;
                if (console) { console->debug("Exit is_distance_less, i = {}, j = {}, spared ub", i, j); }
                return true;
            }
        }

        Real lb = get_cached_lower_bound(i, j);


        if (lb >= value or (strict and lb > value)) {
            m_spared_hera_calls++;
            if (console) { console->debug("Exit is_distance_less, i = {}, j = {}, spared lb", i, j); }
            return false;
        } else {
            lb = lower_bound(i, j);

            if (m_use_geom_lower_bound) {
                Real glb = get_geom_lower_bound(i, j);
                if (glb > lb and glb >= value and lb < value) {
                    m_geom_lower_bound_useful++;
                    lb = glb;
                }
            }

            if (lb > value or (not strict and lb >= value)) {
                m_spared_hera_calls++;
                if (console) { console->debug("Exit is_distance_less, i = {}, j = {}, spared lb", i, j); }
                return false;
            }
        }

        d = get_distance(i, j);
        if (console) { console->debug("Exit is_distance_less, i = {}, j = {}, hera called", i, j); }
        return d < value or (not strict and d <= value);
    }

    bool is_distance_greater(VertexDescriptor i, VertexDescriptor j,
                             Real value, bool strict)
    {
        return not is_distance_less(i, j, value, not strict);
    }


    Real lower_bound(VertexDescriptor i, VertexDescriptor j)
    {
        if (i == j) {
            return 0.0;
        }

        if (console) { console->debug("Enter lower_bound, i = {}, j = {}", i, j); }
        double result = 0.0;

        if (is_distance_cached(i, j, result)) {
            if (console) {
                console->debug("Exit lower_bound, i = {}, j = {}, cached distance, result = {}", i, j, result);
            }
            return result;
        } else {
            return 0.0;
            result = 0.0;
        }

        EdgeIterator ei, ei_end;
        WeightMap weight_map = boost::get(boost::edge_weight, m_graph);

        for (boost::tie(ei, ei_end) = edges(m_graph); ei != ei_end; ++ei) {
            VertexDescriptor u = source(*ei, m_graph);
            VertexDescriptor v = target(*ei, m_graph);
            if (same_component(u, i, m_conn_comps) and same_component(v, j, m_conn_comps)) {
                if (console) {
                    console->debug("In lower_bound, i = {}, j = {}, considering edge u = {}, v = {}", i, j, u, v);
                }
                double ui_dist = upper_bound(u, i);
                double vj_dist = upper_bound(v, j);
                double uj_dist = upper_bound(u, j);
                double vi_dist = upper_bound(v, i);
                double uv_dist = get(weight_map, *ei);
                result = std::max(result, uv_dist - ui_dist - vj_dist);
                result = std::max(result, uv_dist - uj_dist - vi_dist);
                if (console) {
                    console->debug("In lower_bound, i = {}, j = {}, considered edge u = {}, v = {}, result = {}", i, j,
                                   u, v, result);
                }
            }
        }

        if (console) { console->debug("Exit lower_bound, i = {}, j = {}, result = {}", i, j, result); }

        m_lower_bound_cache[std::make_pair(i, j)] = result;
        m_lower_bound_cache[std::make_pair(j, i)] = result;

        return result;
    }

    // 0 : False, 1 : True, -1 : NULL
    int is_distance_less_wo_hera(VertexDescriptor i, VertexDescriptor j, Real value, bool strict)
    {

        if (i == j) {
            if (strict) {
                return (int) (0.0 < value);
            } else {
                return (int) (0.0 <= value);
            }
        }

        m_requested_distances.emplace(i, j);
        m_requested_distances.emplace(j, i);

        int result = -1;
        Real cached_lb = get_cached_lower_bound(i, j);
        Real cached_ub = get_cached_upper_bound(i, j);

        if (cached_lb > value or (not strict and cached_lb >= value)) {
            result = 0;
        } else if (cached_ub < value or (not strict and cached_ub <= value)) {
            result = 1;
        }

        if (-1 == result) {
            Real ub = upper_bound(i, j);
            if (ub < value or (not strict and ub <= value)) {
                result = 1;
            } else {
                Real lb = lower_bound(i, j);
                if (m_use_geom_lower_bound) {
                    Real glb = get_geom_lower_bound(i, j);
                    if (glb > value and lb < value) {
                        m_geom_lower_bound_useful++;
                        lb = glb;
                    }
                }
                if (lb > value or (not strict and lb >= value)) {
                    result = 0;
                }
            }
        }

        if (result >= 0) {
            m_spared_hera_calls++;
        }

        return result;
    }

    int is_distance_greater_wo_hera(VertexDescriptor i, VertexDescriptor j, Real value, bool strict)
    {
        int rev_result = is_distance_less_wo_hera(i, j, value, not strict);
        switch (rev_result) {
            case -1 :
                return -1;
            case 0  :
                return 1;
            case 1  :
                return 0;
        }
    }

    void reset_stats()
    {
        m_hera_distance_calls = 0;
        m_spared_hera_calls = 0;
        m_cached_calls = 0;
        m_geom_lower_bound_useful = 0;
    }

    void clear_cache()
    {
        m_distance_cache.clear();
    }

    std::vector<std::pair<double, int>> get_bounds_quality(const std::vector<double>& deltas)
    {
        std::vector<std::pair<double, int>> result;
        for(double delta : deltas) {
            result.emplace_back(delta, 0);
        }
        for(int i = 0; i < m_num_points; ++i)
            for(int j = i + 1; j < m_num_points; ++j) {
                if (not is_distance_cached(i, j)) {
                    Real lb = lower_bound(i, j);
                    if (lb > 0.0) {
                        Real ub = upper_bound(i, j);
                        if (ub < std::numeric_limits<Real>::infinity()) {
                            Real rel_error = (ub - lb) / lb;
                            for (auto& x : result) {
                                if (rel_error <= x.first)
                                    x.second++;
                            }
                        }
                    }
                }
            }
        return result;
    }



private:

    void cache_distance(VertexDescriptor i, VertexDescriptor j, Real& d)
    {
        Edge e1 = std::make_pair(i, j);
        Edge e2 = std::make_pair(j, i);
        m_lower_bound_cache[e1] = d;
        m_lower_bound_cache[e2] = d;
        m_upper_bound_cache[e1] = d;
        m_upper_bound_cache[e2] = d;
        m_distance_cache[e1] = d;
        m_distance_cache[e2] = d;
    }


    void add_edge(VertexDescriptor i, VertexDescriptor j, Real distance)
    {

        if (console) { console->debug("Enter add_edge, i = {}, j = {}, d = {}", i, j, distance); }
        boost::add_edge(i, j, EdgeWeightProperty(distance), m_graph);
        m_conn_comps.union_set(i, j);
        cache_distance(i, j, distance);
        if (console) { console->debug("Exit add_edge, i = {}, j = {}, d = {}", i, j, distance); }
    }

    bool is_distance_cached(VertexDescriptor i, VertexDescriptor j, Real& d) const
    {
        auto find_res = m_distance_cache.find(std::make_pair(i, j));
        if (find_res != m_distance_cache.end()) {
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


};


template<typename R>
std::ostream& operator<<(std::ostream& os, const typename DynamicSpanner<R>::WassersteinSpacePoint& p)
{
    os << "Point(id = " << p.id << ", dgm.size = " << p.dgm.size() << ") ";
    return os;
}

//template <typename R>
//typename WassersteinSpacePoint<R>::CacheType WassersteinSpacePoint<R>::dist_cache;

//template <typename R>
//size_t WassersteinSpacePoint<R>::n_cache_misses;

#endif // WASSERSTEIN_SPACE_POINT_H

