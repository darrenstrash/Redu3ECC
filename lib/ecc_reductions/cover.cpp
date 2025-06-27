#include <cassert>
#include <utility>
#include <vector>
#include <cstddef>
#include <iostream>
#include <fstream>
#include <algorithm> // for sort

using std::size_t;
using namespace std;

// for hashing
#include <boost/functional/hash.hpp>

#include "cover.hpp"
#include "graph.hpp"
#include "adjacency_list.hpp"

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
    bool is_new_clique = false;
    std::vector<node_t> to_add;
    for (auto it1 = clique.cbegin(); it1 != clique.cend(); it1++) {
        for (auto it2 = std::next(it1, 1); it2 != clique.cend(); it2++) {
            if (not is_covered(*it1, *it2)) {
                is_new_clique = true;
                cover_edge(*it1, *it2);
            }
            if (shadows.contains(*it1)) {
                node_container_t shadow_set = shadows.find(*it1)->second;
                to_add.insert(to_add.end(), shadow_set.begin(), shadow_set.end());
            }
        }
    }

    assert(is_new_clique);

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

bool Cover::verify_cover(ECCGraph const & graph, bool const verbose) const {
    bool verified = true;
    // for each edge in graph, check that it is covered
    edgeset_t all_edges;
    for (auto u : graph.vertices) {
        for (auto v : graph.neighbors(u)) {
            edge_t edge = make_pair(u,v);
            all_edges.insert(edge);
            if (covered_edges.find(edge) == covered_edges.end()) {
                if (verbose) {
                    cout << "VERIFY: FAILED" << endl;
                    cout << "VERIFY: Did not cover edge (" << graph.to_original_id[u] << "," << graph.to_original_id[v] << ")" << endl;
                    verified = false;
                    cout << "NOTE: more errors possible, exiting after first failed check." << endl;
                }
                return false;
            }
        }
    }
    for (auto edge : covered_edges) {
        if (all_edges.find(edge) == all_edges.end()) {
            if (verbose) {
                    cout << "VERIFY: FAILED" << endl;
                    cout << "VERIFY: Consistency issue detected" << endl;
                    cout << "VERIFY: Non-edge (" << graph.to_original_id[get<0>(edge)] << "," << graph.to_original_id[get<1>(edge)] << ") is marked as covered." << endl;
                    verified = false;
                    cout << "NOTE: more errors possible, exiting after first failed check." << endl;
                    return false;
            }
        }
    }

    edgeset_t covered_edges2;

    for (auto const & clique : cliques) {
        for (size_t i = 0; i < clique.size(); i++) {
            for (size_t j = i + 1; j < clique.size(); j++) {
                edge_t edge = make_pair(clique[i], clique[j]);
                if (all_edges.find(edge) == all_edges.end()) {
                    if (verbose) {
                        cout << "VERIFY: FAILED" << endl;
                        cout << "VERIFY: Cover contains non-clique";
                        for (auto u : clique) {
                            cout << " " << graph.to_original_id[u];
                        }
                        cout << endl;
                        cout << "VERIFY: (" << graph.to_original_id[clique[i]] << "," << graph.to_original_id[clique[j]] << ") is not an edge"  << endl;
                        cout << "NOTE: more errors possible, exiting after first failed check." << endl;
                    }
                    return false;
                }

                if (covered_edges.find(edge) == covered_edges.end()) {
                    if (verbose) {
                        cout << "VERIFY: FAILED" << endl;
                        cout << "VERIFY: Consistency issue detected" << endl;
                        cout << "VERIFY: (" << graph.to_original_id[clique[i]] << "," << graph.to_original_id[clique[j]] << ") is not marked covered"  << endl;
                        cout << "VERIFY: but is in clique";
                        for (auto u : clique) {
                            cout << " " << graph.to_original_id[u];
                        }
                        cout << endl;
                        cout << "NOTE: more errors possible, exiting after first failed check." << endl;
                    }
                    return false;
                }
                covered_edges2.insert(edge);
                covered_edges2.insert(make_pair(get<1>(edge), get<0>(edge)));
            }
        }
    }
    if (covered_edges.size() != covered_edges2.size()) {
        if (verbose) {
            cout << "VERIFY: FAILED" << endl;
            cout << "VERIFY: Consistency issue detected" << endl;
            cout << "VERIFY: " << covered_edges.size() << " edges (and their reverses) are marked as covered" << endl;
            cout << "VERIFY: But " << covered_edges2.size() << " edges (and their reverses) are actually covered by cliques..." << endl;
            cout << "NOTE: more errors possible, exiting after first failed check." << endl;
        }
        return false;
    }

    // for each covered edge / clique, check that its edge existed in the graph.
    return true;
}

void Cover::write_cover(ECCGraph const &graph, string const &filename) {
    std::ofstream f(filename.c_str());
    std::cout << "NOTE: Writing cover to " << filename << " ... " << std::endl;

    clique_t original_clique;
    for (auto const & clique : cliques) {
        assert(!clique.empty());
        original_clique.clear();
        original_clique.reserve(clique.size());
        for (node_t v : clique) {
            original_clique.push_back(graph.to_original_id[v]);
        }
        sort(original_clique.begin(), original_clique.end());
        f << original_clique[0];
        for (size_t i = 1; i < original_clique.size(); i++) {
            f << " " << original_clique[i];
        }
        f << endl;
    }
    f.close();
}
