#pragma once

#include <vector>
#include <utility>
#include <string>

#include "cover.hpp"
#include "graph.hpp"

bool is_edge_clique_cover(ECCGraph const& graph, Cover const& cover);
bool decompress_verify(ECCGraph const& graph, Cover const& cover);

void apply_reductions(ECCGraph const& graph, Cover& cover, bool one_enabled, bool two_enabled, bool three_enabled, bool four_enabled, size_t const component = 0);
bool compute_edge_clique_cover(ECCGraph const& graph, Cover& cover, std::string const &basename, size_t const k, size_t& total_calls, bool rule_one_enabled, bool rule_two_enabled, bool rule_three_enabled, bool rule_four_enabled);
void order_neighbors_by_degree(ECCGraph& graph);
void order_neighbors_by_neighbor_degree(ECCGraph& graph);
void compute_connected_components(ECCGraph const& graph, Cover& cover);
