#pragma once

#include <utility> // for pair
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <cstddef>
using std::size_t;

#include "graph.hpp"

struct NodePairHash {
    size_t operator()(std::pair<uint32_t, uint32_t> const& p) const;
};

class Cover {
private:
    std::unordered_set<std::pair<node_t, node_t>, NodePairHash> covered_edges;
    std::unordered_map<node_t, node_container_t> shadows;
public:
    std::vector<node_container_t> cliques;
    std::vector<bool> removed_nodes;
    size_t num_removed_nodes;
    std::vector<size_t> components;
    size_t num_components;
    std::vector<std::vector<node_t>> split_vertices;

    Cover();
    Cover(size_t const);
    Cover(Cover const& other);

    void split_vertex(node_t);
    size_t get_component(node_t) const;
    bool is_covered(node_t, node_t) const;
    void cover_edge(node_t, node_t);
    void cover_clique(node_container_t const&);
    void shadow_node(node_t, node_t);
    bool is_removed(node_t) const;
    void remove_node(node_t);
    size_t num_covered_edges() const;
};
