#include <cassert>
#include <sstream>
#include <string>
#include <istream>

#include "graph.hpp"
#include "adjacency_list.hpp"

ECCGraph::ECCGraph() {
    n = 0;
    e = 0;
}

ECCGraph::ECCGraph(std::istream& is) {
    // assert(ifs.is_open());

    n = 0;
    e = 0;

    std::string line;
    while (std::getline(is, line)) {
        if (line.size() == 0 or line[0] == '#' or line[0] == '%') continue;
        assert(isdigit(line[0]));

        bool in_gap = false;
        bool past_gap = false;
        for (char byte : line) {
            assert(isdigit(byte) or isspace(byte) or byte == ',');
            if (in_gap and isdigit(byte)) {
                in_gap = false;
                assert(not past_gap);
                past_gap = true;
            }
            else if (not in_gap and (isspace(byte) or byte == ',')) {
                in_gap = true;
            }
        }

        std::istringstream iss = std::istringstream(line);

        node_t v1;
        iss >> v1;
        node_t v2;
        iss >> v2;
        if (v1 == v2) continue; // We don't deal with loops

        if (not has_node(v1)) {
            add_node(v1);
        }
        if (not has_node(v2)) {
            add_node(v2);
        }
        add_edge(v1, v2);
    }
}

void ECCGraph::add_node(node_t v) {
    // assert(not has_node(v));

    vertices.push_back(v);
    adj_list.add_node(v);
    n++;
}

void ECCGraph::add_edge(node_t v1, node_t v2) {
    // assert(v1 != v2);

    // assert(not has_edge(v1, v2));

    if (not has_node(v1)) add_node(v1);
    if (not has_node(v2)) add_node(v2);

    if (not adj_list.has_edge(v1, v2)) {
        adj_list.add_edge(v1, v2);
        e++;
    }
}

bool ECCGraph::has_edge(node_t v1, node_t v2) const {
    // assert(has_node(v1));
    // assert(has_node(v2));

    bool result = adj_list.has_edge(v1, v2);
    // assert(result == adj_list.has_edge(v2, v1));

    return result;
}

adj_list_data_t const& ECCGraph::get_adj_list() const {
    return adj_list.data;
}

bool ECCGraph::has_node(node_t v) const {
    return adj_list.has_node(v);
}

node_container_t const& ECCGraph::neighbors(node_t v) const {
    // assert(has_node(v));
    return adj_list.neighbors(v);
}

std::ostream& operator<<(std::ostream& ostream, ECCGraph const& graph) {
    for (auto const& [n, neighbors] : graph.get_adj_list()) {
        ostream << n << ":";
        for (auto neighbor : neighbors) {
            ostream << " " << neighbor;
        }
        ostream << "\n";
    }
    
    return ostream;
}
