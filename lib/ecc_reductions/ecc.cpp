#include <iostream>
#include <vector>
#include <utility> // for pair
#include <limits>
#include <unordered_set>
#include <algorithm> // for sort
#include <cstddef> // for size_t
#include <cstdlib>
#include <fstream>
#include <chrono>

using std::size_t;
using std::vector;

// ECC Reduction / Conversion code
#include "ecc.hpp"
#include "cover.hpp"
#include "graph.hpp"
#include "ecc_to_vcc.hpp"

// KaHIP
#include "graph_io.h"

//bool debug = false;

size_t apply_rule_one(ECCGraph const& graph, Cover& cover, size_t const component = 0) {
    size_t ret = 0;

    size_t vertex_count = 0;

    for (auto const& [v1, v1_neighbors] : graph.get_adj_list()) {
        if (cover.is_removed(v1) or (component != 0 and cover.components[v1] != component)) continue;
        vertex_count++;

        bool all_edges_covered = true;
        for (node_t v2 : v1_neighbors) {
            if (cover.is_removed(v2)) continue;

            if (not cover.is_covered(v1, v2)) {
                all_edges_covered = false;
                break;
            }
        }

        if (all_edges_covered) {
            cover.remove_node(v1);
            //std::cout << "Rule 1 is removing " << v1 << "\n";
            for (node_t neighbor : graph.neighbors(v1)) {
                cover.cover_edge(neighbor, v1);
            }
            ret++;
        }
////        if (vertex_count % 10000 == 0) {
////            std::cerr << "(R1) Edges: " << graph.e << " -> ";
////            std::cerr << graph.e - cover.num_covered_edges() << std::endl;
////        }
    }

    return ret;
}

void compute_common_neighbors(ECCGraph const& graph, Cover const& cover, node_t v1, node_t v2, node_container_t& output) {
    // Assumes nodes contains only valid nodes (not removed)

    for (node_t neighbor : graph.neighbors(v1)) {
        if (cover.is_removed(neighbor) or neighbor == v2) continue;
        if (graph.has_edge(v2, neighbor)) {
            output.insert(neighbor);
        }
    }
}

bool is_clique(ECCGraph const& graph, Cover const& cover, node_container_t const& nodes) {
    for (auto it1 = nodes.cbegin(); it1 != nodes.cend(); it1++) {
        if (cover.is_removed(*it1)) continue;
        for (auto it2 = std::next(it1, 1); it2 != nodes.cend(); it2++) {
            if (cover.is_removed(*it2)) continue;
            if (not graph.has_edge(*it1, *it2)) return false;
        }
    }

    return true;
}

size_t apply_rule_two(ECCGraph const& graph, Cover& cover, size_t const component = 0) {
    size_t ret = 0;
    size_t vertex_count = 0;

    for (auto const& [v1, v1_neighbors] : graph.get_adj_list()) {
        if (cover.is_removed(v1) or (component != 0 and cover.get_component(v1) != component)) continue;
        vertex_count++;
        for (node_t v2 : v1_neighbors) {
            if (v2 <= v1 or cover.is_removed(v2) or cover.is_covered(v1, v2)) continue;

            node_container_t common_neighbors;
            compute_common_neighbors(graph, cover, v1, v2, common_neighbors);
            common_neighbors.insert(v1);
            common_neighbors.insert(v2);

            //debug = (v1 == 790094 and v2 == 53775) or (v1 == 53775 and v2 == 790094);

////            if (debug) {
////                std::cout << "common neighbors of 15 and 45968 = ";
////                for (auto v : common_neighbors) {
////                    std::cout << " " << v;
////                }
////                std::cout << std::endl;
////            }

            if (is_clique(graph, cover, common_neighbors)) {
////                if (debug) {
////                    std::cout << "    is a clique!" << std::endl;
////                }
                cover.cover_clique(common_neighbors);
                // std::cout << "Rule 2 is adding the clique ";
                // for (auto const& n : common_neighbors) std::cout << n << " ";
                // std::cout << "based on edge " << v1 << " " << v2 << "\n";
                ret++;
            } else {
////                if (debug) {
////                    std::cout << "    is not a clique!" << std::endl;
////                }
            }
////            if (vertex_count % 100000 == 0) {
////                std::cerr << "(R2)Edges: " << graph.e << " -> ";
////                std::cerr << graph.e - cover.num_covered_edges() << std::endl;
////            }
        }
    }

    return ret;
}

void compute_prisoners_and_exits(ECCGraph const& graph, Cover const& cover, node_t v1, node_container_t& prisoners, node_container_t& exits) {
    for (node_t v2 : graph.neighbors(v1)) {
        if (cover.is_removed(v2)) continue;

        for (node_t v2_neighbor : graph.neighbors(v2)) {
            if (v2_neighbor == v1 or cover.is_removed(v2_neighbor)) continue;

            if (not graph.has_edge(v1, v2_neighbor)) {
                exits.insert(v2);
                break;
            }
        }
        if (not exits.contains(v2)) {
            prisoners.insert(v2);
        }
    }
}

bool prisoners_dominate_exits(ECCGraph const& graph, Cover const& cover, node_container_t const& prisoners, node_container_t const& exits) {
    // Assumes prisoners and exits all exist.
    if (prisoners.empty()) return false;
    if (exits.empty()) return true;
    
    for (node_t exit : exits) {
        bool dominated = false;
        for (node_t prisoner : prisoners) { // It might be better to loop over the neighbors of the exit.
            if (graph.has_edge(exit, prisoner) and not cover.is_covered(exit, prisoner)) { // This is beyond what the paper said. For justification, try on graphs/easy/right.graph
                dominated = true;
                break;
            }
        }
        if (not dominated) {
            return false;
        }
    }

    return true;
}

size_t apply_rule_three(ECCGraph const& graph, Cover& cover, size_t const component = 0) {

    for (auto const& v : graph.vertices) { // Remove graph.vertices at some point
        if (cover.is_removed(v)) continue;

        node_container_t prisoners;
        node_container_t exits;
        compute_prisoners_and_exits(graph, cover, v, prisoners, exits);

        // Not necessary when fully reduced wrt rules 1 and 2
        bool found_valid_nbr = false;
        for (node_t const& prisoner : prisoners) {
            found_valid_nbr = false;
            for (auto const& nbr_of_prisoner : graph.neighbors(prisoner)) {
                if (nbr_of_prisoner != v and not cover.is_removed(nbr_of_prisoner) and not cover.is_covered(nbr_of_prisoner, prisoner)) {
                    found_valid_nbr = true;
                    break;
                }
            }
            if (not found_valid_nbr) break;
        }
        if (not found_valid_nbr) {
            continue;
        }
        

        if (prisoners_dominate_exits(graph, cover, prisoners, exits)) {
            // std::cout << "Applying rule 3 to " << v << "\n";
            // std::cout << "Prisoners: ";
            // for (auto const& p: prisoners) std::cout << p << " ";
            // std::cout << "\n";
            // std::cout << "Exits: ";
            // for (auto const& e : exits) std::cout << e << " ";
            // std::cout << "\n";

            for (node_t prisoner : prisoners) {
                cover.shadow_node(prisoner, v);
            }
            cover.remove_node(v);
            for (node_t neighbor : graph.neighbors(v)) {
                cover.cover_edge(v, neighbor);
            }
            return true;
        }
    }

    return false;
}

bool is_edge_clique_cover(ECCGraph const& graph, Cover const& cover) {
    return cover.num_covered_edges() == graph.e;
}

bool decompress_verify(ECCGraph const& original, Cover const& cover) {
    ECCGraph decompressed;

    for (node_container_t const& clique: cover.cliques) {
        for (node_t v1 : clique) {
            if (not decompressed.has_node(v1)) {
                decompressed.add_node(v1);
            }

            for (node_t v2 : clique) {
                if (v1 == v2) continue;
                if (not decompressed.has_node(v2)) {
                    decompressed.add_node(v2);
                }
                decompressed.add_edge(v1, v2);
            }
        }
    }

    for (auto const& [v, v_neighbors] : original.get_adj_list()) {
        if (not decompressed.has_node(v)) {
            return false;
        }
        std::unordered_set<node_t> neighbors;
        neighbors.insert(v_neighbors.begin(), v_neighbors.end());
        neighbors.insert(decompressed.neighbors(v).begin(), decompressed.neighbors(v).end());
        if (neighbors.size() != v_neighbors.size() or neighbors.size() != decompressed.neighbors(v).size()) {
            return false;
        }
    }

    return true;
}

void compute_connected_components(ECCGraph const& graph, Cover& cover) {
    // Assumes cover.connected_components is all 0s.
    cover.num_components = 0;
    cover.components = std::vector<size_t>(graph.n, 0);
    size_t seen = 0;
    while (seen != graph.n - cover.num_removed_nodes) {
        std::vector<node_t> static stack;
        stack.clear();
        for (auto const& node : graph.vertices) {
            if (not cover.is_removed(node) and cover.components[node] == 0) {
                stack.push_back(node);
                break;
            }
        }
        cover.num_components++;

        cover.components[stack[stack.size() - 1]] = cover.num_components;
        seen++;
        while (not stack.empty()) {
            node_t const n = stack[stack.size() - 1];
            stack.pop_back();

            for (node_t neighbor : graph.neighbors(n)) {
                if (cover.components[neighbor] == 0) {
                    cover.components[neighbor] = cover.num_components;
                    stack.push_back(neighbor);
                    seen++;
                }
            }
        }
    }
}

void find_cut_vertices(ECCGraph const& graph, Cover const& cover, std::vector<node_t>& cut_vertices) {
    // std::vector<node_t> stack;
    // std::vector<bool> static visited = std::vector<bool>(graph.n, false);               // These will require optimization
    // std::vector<size_t> static lowest_neighbor_time = std::vector<size_t>(graph.n, 0);
    // std::vector<size_t> static discovery_time = std::vector<size_t>(graph.n, 0);
    // std::vector<node_t> static parents = std::vector<node_t>(graph.n, 0);
    // size_t time = 0;
    // size_t children = 0;
    // for (node_t n : graph.vertices) {
    //     if (not cover.is_removed(n)) {
    //         stack.push_back(n);
    //         break;
    //     }
    // }
    // while (not stack.empty()) {
    //     node_t n = stack[stack.size() - 1];
    //     stack.pop_back();
    // 
    //     visited[n] = true;
    //     lowest_neighbor_time[n] = time;
    //     discovery_time[n] = time;
    //     time++;
    // 
    //     for (node_t neighbor : graph.neighbors(n)) {
    //         if (cover.is_removed(neighbor)) continue;
    //         if (not visited[neighbor]) {
    //             parent[neighbor] = n;
    //             children++;
    //             stack.push_back(neighbor);
    //         }
    //     }
    // }
}

bool apply_rule_four(ECCGraph const& graph, Cover& cover, size_t const component = 0) {
    return false;
}

void apply_reductions(ECCGraph const& graph, Cover& cover, bool one_enabled, bool two_enabled, bool three_enabled, bool four_enabled, size_t const component) {
    bool need_to_compute_connected_components = true;
    bool first = true;
    while (true) {
        if (!first and one_enabled and apply_rule_one(graph, cover, component)) {
            need_to_compute_connected_components = true;

////            std::cerr << "    Nodes: " << graph.n << " -> ";
////            size_t remaining_nodes = 0;
////            for (auto const& [v1, v1_neighbors] : graph.get_adj_list()) {
////                if (not cover.is_removed(v1))
////                    remaining_nodes++;
////            }
////            std::cerr << remaining_nodes << std::endl;
            std::cerr << "    Edges: " << graph.e << " -> ";
            std::cerr << graph.e - cover.num_covered_edges() << std::endl;
            continue;
        }
        if (two_enabled and apply_rule_two(graph, cover, component)) {
            first = false;
            continue;
        }
        if (three_enabled and apply_rule_three(graph, cover, component)) {
            first = false;
            need_to_compute_connected_components = true;
            continue;
        }
        // if (four_enabled and need_to_compute_connected_components) {
        //     compute_connected_components(graph, cover);
        //     apply_rule_four(graph, cover, component);
        //     need_to_compute_connected_components = false;
        //     continue;
        // }

        // if (need_to_compute_connected_components) {
        //     compute_connected_components(graph, cover);
        // }
        return;
    }
}

void compute_all_maximal_cliques(ECCGraph const& graph, Cover const& cover, node_container_t& R, node_container_t& P, node_container_t& X, std::vector<node_container_t>& cliques) {
    // This is just Bron Kerbosch
    if (P.empty() and X.empty()) {
        cliques.push_back(R);
        // std::cerr << "Got clique " << cliques.size() << "\n";
        return;
    }

    while (not P.empty()) {
        node_t const v = *P.cbegin();
        node_container_t new_R = R;
        new_R.insert(v);
        node_container_t new_P;
        node_container_t new_X;
        for (node_t n : graph.neighbors(v)) {
            if (cover.is_removed(n)) continue;
            if (P.contains(n)) {
                new_P.insert(n);
            }
            if (X.contains(n)) {
                new_X.insert(n);
            }
        }
        compute_all_maximal_cliques(graph, cover, new_R, new_P, new_X, cliques);
        P.erase(v);
        X.insert(v);
    }
}

size_t compute_score(ECCGraph const& graph, Cover const& cover, node_t v1, node_t v2) {
    size_t result = 0;
    node_container_t common_neighbors;
    compute_common_neighbors(graph, cover, v1, v2, common_neighbors);

    for (auto it1 = common_neighbors.cbegin(); it1 != common_neighbors.cend(); it1++) {
        if (cover.is_removed(*it1)) continue;
        for (auto it2 = std::next(it1, 1); it2 != common_neighbors.cend(); it2++) {
            if (cover.is_removed(*it2)) continue;
            result += not graph.has_edge(*it1, *it2);
        }
    }

    return result;
}

std::pair<node_t, node_t> pick_lowest_score_edge(ECCGraph const& graph, Cover const& cover) {
    // Pick the lowest score edge. (The edge whose common neighbors are closest to forming a clique)
    size_t best_score = std::numeric_limits<size_t>::max();
    std::pair<node_t, node_t> best_edge = {std::numeric_limits<node_t>::max(), std::numeric_limits<node_t>::max()};

    for (auto const& [v1, v1_neighbors] : graph.get_adj_list()) {
        if (cover.is_removed(v1)) continue;
        for (node_t v2 : v1_neighbors) {
            if (v2 <= v1 or cover.is_removed(v2) or cover.is_covered(v1, v2)) continue;
            size_t score = compute_score(graph, cover, v1, v2);
            if (score < best_score) {
                best_score = score;
                best_edge = {v1, v2};
            }
        }
    }

    return best_edge;
}

std::pair<node_t, node_t> pick_highest_score_edge(ECCGraph const& graph, Cover const& cover) {
    // Pick the lowest score edge. (The edge whose common neighbors are closest to forming a clique)
    size_t best_score = 0;
    std::pair<node_t, node_t> best_edge = {std::numeric_limits<node_t>::max(), std::numeric_limits<node_t>::max()};

    for (auto const& [v1, v1_neighbors] : graph.get_adj_list()) {
        if (cover.is_removed(v1)) continue;
        for (node_t v2 : v1_neighbors) {
            if (v2 <= v1 or cover.is_removed(v2) or cover.is_covered(v1, v2)) continue;
            size_t score = compute_score(graph, cover, v1, v2);
            if (score > best_score) {
                best_score = score;
                best_edge = {v1, v2};
            }
        }
    }

    return best_edge;
}

std::pair<node_t, node_t> pick_first_edge(ECCGraph const& graph, Cover const& cover) {
    // Pick the first edge.
    for (auto const& [v1, v1_neighbors] : graph.get_adj_list()) {
        if (cover.is_removed(v1)) continue;
        for (node_t v2 : v1_neighbors) {
            if (cover.is_removed(v2) or cover.is_covered(v1, v2)) continue;
            return {v1, v2};
        }
    }
    return {std::numeric_limits<node_t>::max(), std::numeric_limits<node_t>::max()};
}

void order_neighbors_by_neighbor_degree(ECCGraph& graph) {
    // Sort the adjacency lists so the nodes with the smallest 2-neighborhood show up first.
    
    auto keyfn = [graph](node_t const& a, node_t const& b) {
        size_t a_sum = 0;
        size_t b_sum = 0;
        for (auto neighbor_neighbor : graph.neighbors(a)) {
            a_sum += graph.neighbors(neighbor_neighbor).size();
        }
        for (auto neighbor_neighbor : graph.neighbors(b)) {
            b_sum += graph.neighbors(neighbor_neighbor).size();
        }
        return a_sum > b_sum;
    };
    
    // std::sort(graph.vertices.begin(), graph.vertices.end(), keyfn); // Makes things way too slow on big graphs

    for (auto& [n1, neighbors] : graph.adj_list.data) {
        std::sort(neighbors.begin(), neighbors.end(), keyfn);
    }
}

void order_neighbors_by_degree(ECCGraph& graph) {
    // Sort the adjacency lists so the nodes with the smallest 2-neighborhood show up first.
    
    auto keyfn = [graph](node_t const& a, node_t const& b) {
        return graph.neighbors(a).size() > graph.neighbors(b).size();
    };
    
    // std::sort(graph.vertices.begin(), graph.vertices.end(), keyfn); // Makes things way too slow on big graphs
    
    for (auto& [n1, neighbors] : graph.adj_list.data) {
        std::sort(neighbors.begin(), neighbors.end(), keyfn);
    }
}

size_t size_of_big_vertex_independent_set(ECCGraph const& graph, Cover const& cover) {
    // Finds a large independent set among the uncovered vertices in the graph.

    std::vector<node_t> ind_set;
    for (auto const& [n1, _] : graph.get_adj_list()) {
        if (cover.is_removed(n1)) continue;
        bool independent = true;
        for (auto const& n2 : ind_set) {
            // Don't need to check if removed here, because we only add valid nodes to ind_set.
            if (graph.has_edge(n1, n2)) {
                independent = false;
                break;
            }
        }
        if (independent) ind_set.push_back(n1);
    }
    return ind_set.size();
}

size_t size_of_big_edge_independent_set(ECCGraph const& graph, Cover const& cover) {
    // Computes a large set of edges such that the subgraph induced by their vertices is K_4-free.

    std::vector<std::pair<node_t, node_t>> ind_set;
    for (auto const& n1 : graph.vertices) {
        if (cover.is_removed(n1)) continue;
        
        auto const& neighbors = graph.neighbors(n1);
        bool independent = true;

        for (auto const& n2 : neighbors) {
            if (n1 >= n2 or cover.is_removed(n2) or cover.is_covered(n1, n2)) continue;
            for (auto const& [v1, v2] : ind_set) {
                if ((n1 == v1 and graph.has_edge(v2, n2)) or
                    (n2 == v2 and graph.has_edge(v1, n1)) or
                    (n1 == v2 and graph.has_edge(v1, n2)) or
                    (n2 == v1 and graph.has_edge(v2, n1)) or
                    (graph.has_edge(n1, v1) and graph.has_edge(n2, v1) and graph.has_edge(n1, v2) and graph.has_edge(n2, v2))) {
                    independent = false;
                    break;
                }
            }
            if (independent) {
                ind_set.push_back({n1, n2});
                break;
            }
        }
    }
    return ind_set.size();
}

bool compute_edge_clique_cover(ECCGraph const& graph, Cover& cover, std::string const &basename, size_t const k, size_t& total_calls, bool rule_one_enabled, bool rule_two_enabled, bool rule_three_enabled, bool rule_four_enabled) {
    total_calls++;

    std::cerr << "Applying reductions\n";

    auto const start = std::chrono::high_resolution_clock::now();
    apply_reductions(graph, cover, rule_one_enabled, rule_two_enabled, rule_three_enabled, rule_four_enabled);

    auto const end = std::chrono::high_resolution_clock::now();

    std::cerr << "    Nodes: " << graph.n << " -> ";
    size_t remaining_nodes = 0;
    for (auto const& [v1, v1_neighbors] : graph.get_adj_list()) {
        if (not cover.is_removed(v1))
            remaining_nodes++;
    }
    std::cerr << remaining_nodes << std::endl;
    std::cerr << "    Edges: " << graph.e << " -> ";
    std::cerr << graph.e - cover.num_covered_edges() << std::endl;

    std::chrono::duration<double, std::milli> const time_to_solution = end - start;

    //std::cerr << "time_to_solution=" << time_to_solution.count() << std::endl;

/**
    std::string kernel_all = basename + ".kernel.all.edgelist";
    std::string kernel_uncovered = basename + ".kernel.uncovered.edgelist";

    std::ofstream fall(kernel_all);
    std::ofstream funcovered(kernel_uncovered);

    size_t kernel_vertices = 0;
    size_t kernel_all_edges = 0;
    size_t kernel_uncovered_edges = 0;

    for (auto const& [v1, v1_neighbors] : graph.get_adj_list()) {
        if (cover.is_removed(v1)) continue;
        kernel_vertices++;
        for (node_t v2 : v1_neighbors) {
            if (cover.is_removed(v2)) continue;
            if (not cover.is_covered(v1, v2)) {
                funcovered << v1 << " " << v2 << std::endl;
                kernel_uncovered_edges++;
            }
            fall << v1 << " " << v2 << std::endl;
            kernel_all_edges++;
        }
    }

    fall.close();
    funcovered.close();
**/

    std::cerr << "input_vertices=" << graph.n << std::endl;
    std::cerr << "input_edges=" << graph.e / 2 << std::endl;

    ECC2VCC converter(graph, cover);
    vector<vector<NodeID>> vcc_adjlist;
    vcc_adjlist = converter.ecc_to_vcc();

    //std::cerr << "kernel_vertices=" << kernel_vertices << std::endl;
    //std::cerr << "kernel_all_edges=" << kernel_all_edges / 2 << std::endl;
    //std::cerr << "kernel_uncovered_edges=" << kernel_uncovered_edges / 2 << std::endl;
    std::cerr << "num_cliques=" << cover.cliques.size() << std::endl;
    std::cerr << "time_to_solution=" << time_to_solution.count() / 1000.0 << std::endl;

    exit(0);
    if (cover.cliques.size() >= k) return false;
    if (is_edge_clique_cover(graph, cover)) return true;

    // size_t lower_bound = cover.cliques.size() + size_of_big_edge_independent_set(graph, cover);
    // std::cerr << "Found a lower bound of " << lower_bound << "\n";
    // if (lower_bound >= k) return false;

    // Pick edge to branch on
    std::pair<node_t, node_t> const lowest_score_edge = pick_lowest_score_edge(graph, cover);
    // std::cerr << "Got lowest score edge.\n";
    if (lowest_score_edge.first == std::numeric_limits<node_t>::max()) {
        std::cerr << "Best edge never got set. Probably every vertex is deleted.\n";
        size_t removed_vertices = 0;
        for (bool removed : cover.removed_nodes) removed_vertices += removed;
        std::cerr << removed_vertices << "/" << graph.n << " vertices have been removed, and " << cover.num_covered_edges() << "/" << graph.e << " edges have been covered.\n";
        exit(1);
    }

    // R contains vertices that have to be in the clique (for us, the vertices in the edge we're branching on)
    node_container_t static R;
    R.clear();
    R.insert(lowest_score_edge.first);
    R.insert(lowest_score_edge.second);

    // Contains vertices that could be in the clique (for us, the common neighbors)
    node_container_t static common_neighbors;
    common_neighbors.clear();
    compute_common_neighbors(graph, cover, lowest_score_edge.first, lowest_score_edge.second, common_neighbors);
    // std::cerr << "Got common neighbors.\n";
    // X contains the clique as it is built. It starts empty.
    node_container_t static X;
    X.clear();
    std::vector<node_container_t> lowest_score_edge_cliques; // Cannot be static because used in multiple recursive calls at once.
    compute_all_maximal_cliques(graph, cover, R, common_neighbors, X, lowest_score_edge_cliques);
    
    Cover best_cover = Cover(graph.n);
    bool found_better_cover = false;
    size_t new_k = k;

    std::vector<node_container_t>& cliques = lowest_score_edge_cliques; // Left for future experiments with edge selection

    // std::cout << "Branching on " << lowest_score_edge.first << " " << lowest_score_edge.second << "\n";
    for (node_container_t const& clique : cliques) {
        Cover new_cover = cover; // a copy, on purpose
        new_cover.cover_clique(clique);
        // std::cerr << "Starting branch with this clique: ";
        // for (auto const& n : clique) std::cerr << n << " ";
        // std::cerr << "\n";
        
        bool result = compute_edge_clique_cover(graph, new_cover, basename, new_k, total_calls, rule_one_enabled, rule_two_enabled, rule_three_enabled, rule_four_enabled);
        if (result and new_cover.cliques.size() < new_k) {
            new_k = new_cover.cliques.size();
            best_cover = std::move(new_cover);
            found_better_cover = true;
            // std::cerr << "New best cover is " << new_cover.cliques.size() << "\n";
        }
        // std::cout << "Branch with clique";
        // for (auto const& n : clique) std::cout << n << " ";
        // std::cout << "over.";
    }

    if (found_better_cover) {
        cover = std::move(best_cover);
    }

    return is_edge_clique_cover(graph, cover);
}
