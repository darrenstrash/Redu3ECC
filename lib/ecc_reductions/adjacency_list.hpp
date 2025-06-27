#pragma once

#include "types.hpp"

class AdjacencyList {
public:
    adj_list_data_t data; // Ideally this wouldn't be public. I'll fix this later.

    AdjacencyList();

    void add_edge(node_t, node_t);
    void add_node(node_t);
    bool has_edge(node_t, node_t) const;
    bool has_node(node_t) const;
    node_container_t const& neighbors(node_t) const;

////private:
////    std::unordered_set<std::pair<node_t, node_t>, NodePairHash2> edges;
};
