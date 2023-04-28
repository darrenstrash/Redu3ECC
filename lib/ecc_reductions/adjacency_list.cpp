#include <algorithm>
#include "adjacency_list.hpp"
#include <initializer_list>

// for hashing
#include <boost/functional/hash.hpp>

////size_t NodePairHash2::operator()(std::pair<uint32_t, uint32_t> const& p) const {
////    // Ensure that size_t, the type of the hash, is large enough
////    // assert(sizeof(size_t) >= sizeof(uint32_t) * 2); // It usually is
////    //return (((size_t)p.first) << sizeof(uint32_t)) | (size_t)p.second;
////    //DS: a real hash function
////    size_t seed = 0;
////    boost::hash_combine(seed, std::get<0>(p));
////    boost::hash_combine(seed, std::get<1>(p));
////    return seed;
////}
////

// please don't call this
bool NodeVector::contains(node_t const& node) const {
    assert(0);
    return std::find(cbegin(), cend(), node) != cend();
}

void NodeVector::insert(node_t const& node) {
    push_back(node);
}

// please don't call this either
void NodeVector::erase(node_t const& node) {
    assert(0);
    // two full passes. Should just perform one pass, swap to end and remove.
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
////        edges.insert({v1, v2});
////        edges.insert({v2, v1});
    }
}

void AdjacencyList::add_node(node_t v) {
    data.insert({v, {}});
}

bool AdjacencyList::has_edge(node_t v1, node_t v2) const {
    return data.find(v1)->second.contains(v2) or data.find(v2)->second.contains(v1);
////    return edges.count({v1, v2});
}

bool AdjacencyList::has_node(node_t v) const {
    return data.contains(v);
}

node_container_t const& AdjacencyList::neighbors(node_t v) const {
    return data.find(v)->second;
}
