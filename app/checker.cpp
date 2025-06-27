/******************************************************************************
 * checker.cpp 
 * *
 * Source of Redu3ECC
 * Darren Strash
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <limits>
#include <vector>
#include <unordered_set>

#include "ecc_reductions/types.hpp"

using namespace std;

// this program implements the functions to check the metis graph 
// format
int main(int argn, char **argv)
{

    if (argn < 2 || argn > 4) {
        std::cout <<  "Usage: checker GRAPH_FILE [COVER_FILE]"  << std::endl;
        exit(0);
    }

    std::cout <<  "*******************************************************************************"  << std::endl;
    std::cout <<  "Redu3ECC -- Graph and cover checker."  << std::endl;
    std::cout <<  "*******************************************************************************"  << std::endl;

    std::string line;
    std::string graphname(argv[1]);

    std::cout << "Reading graph " << graphname << std::endl;
    // open file for reading
    std::ifstream in(graphname.c_str());
    if (!in) {
        std::cerr << "Error opening " << graphname << std::endl;
        return 1;
    }


    edgeset_t edges;
    size_t num_edges = 0;
    while (std::getline(in, line)) {
        if (line.size() == 0 or line[0] == '#' or line[0] == '%') continue;
        assert(isdigit(line[0]));

        bool in_gap = false;
        bool past_gap = false;
        for (char byte : line) {
            if (!isdigit(byte) and !isspace(byte)) {
                std::cout << "ERROR: expecting a line to consist of non-negative integers separated by spaces, got \"" << line << "\"" << std::endl;
                std::cout << "NOTE: More errors possible, exiting after first failed check." << std::endl;
                return 1;
            }
            if (in_gap and isdigit(byte)) {
                in_gap = false;
                if (past_gap) {
                    std::cout << "ERROR: expecting a line to consist of a pair of non-negative integers separated by spaces, got \"" << line << "\"" << std::endl;
                    std::cout << "NOTE: More errors possible, exiting after first failed check." << std::endl;
                    return 1;
                }
                past_gap = true;
            }
            else if (not in_gap and isspace(byte)) {
                in_gap = true;
            }
        }

        std::istringstream iss = std::istringstream(line);

        size_t v1;
        iss >> v1;
        size_t v2;
        iss >> v2;

        if (v1 > numeric_limits<node_t>::max()) {
            std::cout << "ERROR: input contains vertex that does not fit in 32 bits: " << v1 << std::endl;
            std::cout << "NOTE: More errors possible, exiting after first failed check." << std::endl;
            return 1;
        }

        if (v2 > numeric_limits<node_t>::max()) {
            std::cout << "ERROR: input contains vertex that does not fit in 32 bits: " << v2 << std::endl;
            std::cout << "NOTE: More errors possible, exiting after first failed check." << std::endl;
            return 1;
        }

        if (v1 == v2) {
            std::cout << "WARNING: Graph contains loop (" << v1 << "," << v1 << "). Removing..." << std::endl;
            continue; // We don't deal with loops
        }

        edge_t edge =  make_pair(static_cast<node_t>(v1),static_cast<node_t>(v2));
        if (edges.find(edge) != edges.end()) {
            std::cout << "WARNING: Graph contains parallel edge (" << v1 << "," << v1 << "). Removing..." << std::endl;
            continue; // We don't deal with parallel edges
        }

        num_edges++;

        edges.insert(edge);
    }

    std::cout <<  "Graph IO done. Checking graph symmetry"  << std::endl;

    for (edge_t const & edge : edges) {
        edge_t reverse_edge = make_pair(get<1>(edge), get<0>(edge));
        if (edges.find(reverse_edge) == edges.end()) {
            std::cout << "ERROR: Graph is asymmetric:" << std::endl;
            std::cout << "ERROR:         edge (" << get<0>(edge) << "," << get<1>(edge) << ") present, but" << std::endl;
            std::cout << "ERROR: reverse edge (" << get<0>(reverse_edge) << "," << get<1>(reverse_edge) << ") is not." << std::endl;
            std::cout << "NOTE: More errors possible, exiting after first failed check." << std::endl;
            return 1;
        }
    }

    std::cout <<  "The graph format seems correct."  << std::endl;

    if (argn < 3) {
        return 0;
    }

    std::string covername(argv[2]);
    std::cout << "Reading cover file " << covername << std::endl;
    // open file for reading
    std::ifstream in2(covername.c_str());
    if (!in2) {
        std::cerr << "Error opening " << covername << std::endl;
        return 1;
    }

    edgeset_t covered_edges;
    size_t num_cliques = 0;
    while (std::getline(in2, line)) {
        if (line.size() == 0 or line[0] == '#' or line[0] == '%') continue;
        assert(isdigit(line[0]));

        for (char byte : line) {
            if (!isdigit(byte) and !isspace(byte)) {
                std::cout << "ERROR: expecting a line to consist of non-negative integers separated by spaces, got \"" << line << "\"" << std::endl;
                std::cout << "NOTE: More errors possible, exiting after first failed check." << std::endl;
                return 1;
            }
        }

        num_cliques++;
        std::istringstream iss = std::istringstream(line);

        clique_t clique;
        size_t v;
        while (iss.rdstate() != std::ifstream::eofbit) {
            iss >> v;
            if (v > numeric_limits<node_t>::max()) {
                std::cout << "ERROR: input contains vertex that does not fit in 32 bits: " << v << std::endl;
                std::cout << "NOTE: More errors possible, exiting after first failed check." << std::endl;
                return 1;
            }
            clique.push_back(static_cast<node_t>(v));
        }

        clique_t sorted_clique = clique;
        sort(sorted_clique.begin(), sorted_clique.end());
        for (size_t i = 1; i < sorted_clique.size(); i++) {
            if (sorted_clique[i-1] == sorted_clique[i]) {
                std::cout << "WARNING: Clique ";
                for (node_t v : clique) {
                    std::cout << " " << v;
                }
                std::cout << " contains duplicate vertex " << sorted_clique[i-1] << ". Removing..." << std::endl;
            }
        }

        auto it = std::unique(sorted_clique.begin(), sorted_clique.end());
        sorted_clique.resize(std::distance(sorted_clique.begin(), it));
        //cover.emplace_back(sorted_clique);

        if (sorted_clique.size() == 1) {
            std::cout << "WARNING: Clique ";
            for (node_t v : clique) {
                std::cout << " " << v;
            }
            std::cout << " covers no edges. Removing from consideration..." << std::endl;
            continue;
        }

        for (size_t i = 0; i < sorted_clique.size(); i++) {
            for (size_t j = i + 1; j < sorted_clique.size(); j++) {
                edge_t edge = make_pair(sorted_clique[i], sorted_clique[j]);
                if (edges.find(edge) == edges.end()) {
                    std::cout << "ERROR: ";
                    for (node_t v : clique) {
                        std::cout << " " << v;
                    }
                    std::cout << " is not a clique. Contains non-edge (" << get<0>(edge) << "," << get<1>(edge) << ")." << std::endl;
                    std::cout << "NOTE: More errors possible, exiting after first failed check." << std::endl;
                    return 1;
                }
                covered_edges.insert(edge);
                covered_edges.insert(make_pair(get<1>(edge), get<0>(edge)));
            }
        }
    }

    if (edges.size() != covered_edges.size()) {
        std::cout << "ERROR: At least one edge is uncovered: ";
        for (edge_t edge : edges) {
            if (covered_edges.find(edge) == covered_edges.end()) {
                std::cout << "(" << get<0>(edge) << "," << get<1>(edge) << ")" << std::endl;
                break;
            }
        }
        std::cout << "NOTE: More errors possible, exiting after first failed check." << std::endl;
        return 1;
    }

    std::cout << "The cover has " << num_cliques << " cliques, and covers all " << edges.size() / 2 <<" edges." << std::endl;
    return 0;
}

