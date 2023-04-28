#pragma once

#include <unordered_set>
#include <vector>
#include <string>
#include <istream>

#include "adjacency_list.hpp"

class ECCGraph {
public:
    AdjacencyList adj_list;
    size_t n;
    size_t e;
    std::vector<node_t> vertices;

    ECCGraph();
    ECCGraph(std::istream&);

    void write_edge_list_file(std::string const&) const;

    adj_list_data_t const& get_adj_list() const;

    // These methods do precondition checks, then call corresponding methods in AdjacencyList.
    void add_edge(node_t, node_t);
    void add_node(node_t);
    bool has_edge(node_t, node_t) const;
    bool has_node(node_t v) const;
    node_container_t const& neighbors(node_t) const;
};

std::ostream& operator<<(std::ostream& ostream, ECCGraph const& graph);
