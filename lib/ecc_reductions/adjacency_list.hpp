#pragma once

#include <unordered_map>
#include <cstdint>
#include <vector>

typedef uint32_t node_t; // I might eventually need more than 4 billion nodes...

class NodeVector : public std::vector<node_t> {
public:
    bool contains(node_t const&) const;
    void insert(node_t const&);
    void erase(node_t const&);
};

typedef NodeVector node_container_t;
typedef std::unordered_map<node_t, node_container_t> adj_list_data_t;

class AdjacencyList {
public:
    adj_list_data_t data; // Ideally this wouldn't be public. I'll fix this later.

    AdjacencyList();

    void add_edge(node_t, node_t);
    void add_node(node_t);
    bool has_edge(node_t, node_t) const;
    bool has_node(node_t) const;
    node_container_t const& neighbors(node_t) const;
};