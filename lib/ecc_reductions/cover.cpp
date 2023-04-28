#include <cassert>
#include <utility>
#include <vector>
#include <cstddef>
using std::size_t;

#include "cover.hpp"
#include "graph.hpp"
#include "adjacency_list.hpp"

size_t NodePairHash::operator()(std::pair<uint32_t, uint32_t> const& p) const {
    // Ensure that size_t, the type of the hash, is large enough
    // assert(sizeof(size_t) >= sizeof(uint32_t) * 2); // It usually is
    return (((size_t)p.first) << sizeof(uint32_t)) | (size_t)p.second;
}

Cover::Cover() {
    num_components = 0;
    num_removed_nodes = 0;
}

Cover::Cover(size_t const number_of_nodes) {
    num_components = 0;
    num_removed_nodes = 0;
    removed_nodes = std::vector<bool>(number_of_nodes, false);
    components = std::vector<size_t>(number_of_nodes, 0);
}

Cover::Cover(Cover const& other) {
    num_removed_nodes = other.num_removed_nodes;
    num_components = other.num_components;
    removed_nodes = other.removed_nodes;
    components = other.components;
    covered_edges = other.covered_edges;
    shadows = other.shadows;
    cliques = other.cliques;
}

size_t Cover::get_component(node_t n) const {
    return components[n];
}

void Cover::split_vertex(node_t n) {
    // Requires that connected components have already been computed.
    if (not removed_nodes[n]) {
        removed_nodes[n] = true;
        split_vertices.push_back({n});
    }
}

bool Cover::is_covered(node_t v1, node_t v2) const {
    return covered_edges.contains({v1, v2});
}

void Cover::remove_node(node_t v) {
    if (not removed_nodes[v]) {
        num_removed_nodes++;
    }
    removed_nodes[v] = true;
}

void Cover::cover_edge(node_t v1, node_t v2) {
    if (not covered_edges.contains({v1, v2})) {
        covered_edges.insert({v1, v2});
    }
    if (not covered_edges.contains({v2, v1})) {
        covered_edges.insert({v2,v1});
    }
}

void Cover::shadow_node(node_t target, node_t shadow) {
    // Shadow will be inserted into any clique that target is in.
    if (not shadows.contains(target)) {
        shadows.insert({target, {}});
    }
    shadows.find(target)->second.insert(shadow);
}

void Cover::cover_clique(node_container_t const& clique) {
    // bool is_new_clique = false;
    std::vector<node_t> to_add;
    for (auto it1 = clique.cbegin(); it1 != clique.cend(); it1++) {
        for (auto it2 = std::next(it1, 1); it2 != clique.cend(); it2++) {
            if (not is_covered(*it1, *it2)) {
                // is_new_clique = true;
                cover_edge(*it1, *it2);
            }
            if (shadows.contains(*it1)) {
                node_container_t shadow_set = shadows.find(*it1)->second;
                to_add.insert(to_add.end(), shadow_set.begin(), shadow_set.end());
            }
        }
    }

    for (size_t i = 0; i < to_add.size(); i++) {
        if (shadows.contains(to_add[i])) { // If the things we just added have shadows,
            node_container_t shadow_set = shadows.find(to_add[i])->second;
            to_add.insert(to_add.end(), shadow_set.begin(), shadow_set.end()); // don't forget about them.
                                                                               // (These also get handled in this loop)
                                                                               // I think circular shadows are impossible, so that shouldn't be an issue.
        }
    }

    // assert(is_new_clique);
    cliques.push_back(std::move(clique));
    for (auto const& node : to_add) {
        cliques[cliques.size() - 1].insert(node);
    }
}

bool Cover::is_removed(node_t v) const {
    return removed_nodes[v];
}

size_t Cover::num_covered_edges() const {
    return covered_edges.size() / 2;
}
