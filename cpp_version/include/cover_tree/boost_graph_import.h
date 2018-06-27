#ifndef SPANNER_WASSERSTEIN_BOOST_GRAPH_IMPORT_H
#define SPANNER_WASSERSTEIN_BOOST_GRAPH_IMPORT_H

#include <boost/config.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/incremental_components.hpp>

using boost::graph_traits;
using boost::adjacency_list;
using boost::disjoint_sets;
using boost::property_map;
using boost::edge_weight_t;
using boost::property;
using boost::no_property;
using boost::vecS;
using boost::undirectedS;
using boost::initialize_incremental_components;
using boost::incremental_components;
using boost::dijkstra_shortest_paths;


using EdgeWeightProperty = property<edge_weight_t, double>;
using Graph = adjacency_list<vecS, vecS, undirectedS, no_property, EdgeWeightProperty>;
using VertexDescriptor = typename graph_traits<Graph>::vertex_descriptor;
using EdgeDescriptor = typename graph_traits<Graph>::edge_descriptor;
using EdgeIterator = typename graph_traits<Graph>::edge_iterator;
using Edge = typename std::pair<VertexDescriptor, VertexDescriptor>;
using DisjointSets = boost::disjoint_sets<size_t*, size_t*>;
using WeightMap = boost::property_map<Graph, edge_weight_t>::type;


struct EdgeHashStruct
{
    size_t operator()(const Edge& p) const
    {
        size_t seed = 0;
        boost::hash_combine(seed, p.first);
        boost::hash_combine(seed, p.second);
        return seed;
    }
};


#endif //SPANNER_WASSERSTEIN_BOOST_GRAPH_IMPORT_H
