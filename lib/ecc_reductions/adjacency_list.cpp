#include <algorithm>
#include "adjacency_list.hpp"
#include <initializer_list>

bool NodeVector::contains(node_t const& node) const {
    return std::find(cbegin(), cend(), node) != cend();
}

void NodeVector::insert(node_t const& node) {
    push_back(node);
}

void NodeVector::erase(node_t const& node) {
    auto it = std::find(cbegin(), cend(), node);
    if (it != cend()) {
        std::vector<node_t>& vec = *this;
        vec.erase(it);
    }
}

// All this find->second business needs to get cleaned up.
AdjacencyList::AdjacencyList() {}

void AdjacencyList::add_edge(node_t v1, node_t v2) {
    if (not has_edge(v1, v2)) {
        data.find(v1)->second.insert(v2);
        data.find(v2)->second.insert(v1);
    }
}

void AdjacencyList::add_node(node_t v) {
    data.insert({v, {}});
}

bool AdjacencyList::has_edge(node_t v1, node_t v2) const {
    return data.find(v1)->second.contains(v2) or data.find(v2)->second.contains(v1);
}

bool AdjacencyList::has_node(node_t v) const {
    return data.contains(v);
}

node_container_t const& AdjacencyList::neighbors(node_t v) const {
    return data.find(v)->second;
}