#include <iostream>
#include <utility> // for std::pair
#include <algorithm>

// containers
#include <vector>
#include <unordered_map>
#include <unordered_set>

// local includes
#include "ecc_to_vcc.hpp"
#include "graph.hpp"

namespace std
{
    template<> struct hash<std::pair<NodeID,NodeID>>
    {
        size_t operator()(std::pair<NodeID,NodeID> const & node_pair) const noexcept
        {
            size_t seed = 0;
            boost::hash_combine(seed, std::get<0>(node_pair));
            boost::hash_combine(seed, std::get<1>(node_pair));
            return seed;
        }
    };
};


using namespace std;

vector<vector<NodeID>> ECC2VCC::compute_ecc_adjlist(std::unordered_set<std::pair<NodeID, NodeID>> & uncovered) const {
    size_t kernel_vertices = 0;
    size_t kernel_all_edges = 0;
    size_t kernel_uncovered_edges = 0;

    vector<vector<NodeID>> ecc_adjlist;
    // TODO: Could be really expensive, better way to reserve?
    //ecc_adjlist.reserve(graph.n);

    // TODO: make new_id = 0; doing this for testing and comparing with existing script

    // step 1: reconstruct adjacency list of ecc kernel
    NodeID new_id = 1;
    std::unordered_map<NodeID, NodeID> to_new_id;

    for (auto const& [v1, v1_neighbors] : graph.get_adj_list()) {
        if (cover.is_removed(v1)) continue;
        kernel_vertices++;

        if (to_new_id.find(v1) == to_new_id.end()) {
            if (new_id == 15 or new_id == 45968)
                std::cout << v1 << " -> " << new_id << std::endl;
            to_new_id[v1] = new_id++;
        }

        for (node_t v2 : v1_neighbors) {

            //if (v2 < v1) continue; // avoid double processing
            if (cover.is_removed(v2)) continue;

            // TODO: move below ifs, only for testing
            if (to_new_id.find(v2) == to_new_id.end()) {
                if (new_id == 15 or new_id == 45968)
                    std::cout << v2 << " -> " << new_id << std::endl;
                to_new_id[v2] = new_id++;
            }

            if (to_new_id[v1] >= ecc_adjlist.capacity()) {
                ecc_adjlist.reserve(ecc_adjlist.capacity() * 2);
            }

            if (to_new_id[v1] >= ecc_adjlist.size()) {
                ecc_adjlist.resize(to_new_id[v1] + 1);
            }

            ecc_adjlist[to_new_id[v1]].push_back(to_new_id[v2]);
            //ecc_adjlist[to_new_id[v2]].push_back(to_new_id[v1]);
            if (not cover.is_covered(v1, v2)) {
                //funcovered << v1 << " " << v2 << std::endl;
                uncovered.insert(std::make_pair(to_new_id[v1], to_new_id[v2]));
                //uncovered.insert(std::make_pair(v2, v1));
                kernel_uncovered_edges++;
            }
            //fall << v1 << " " << v2 << std::endl;
            kernel_all_edges++;
        }
    }

    //std::cout << "old ids id[15]=" << to_new_id[15] << ", id[45968]=" << to_new_id[45968] << std::endl;

    std::cerr << "ecc_kernel_vertices=" << kernel_vertices << std::endl;
    std::cerr << "ecc_kernel_all_edges=" << kernel_all_edges / 2 << std::endl;
    std::cerr << "ecc_kernel_uncovered_edges=" << kernel_uncovered_edges / 2 << std::endl;

    /**
    std::cout << kernel_vertices << " " << kernel_all_edges / 2 << std::endl;

    for (auto i : ecc_adjlist) {
        std::sort(i.begin(), i.end());
        if (!i.empty()) {
            std::cout << i[0];
        }
        for (size_t j = 1; j < i.size(); j++) {
            std::cout << " " << i[j];
        }
        std::cout << std::endl;
    }
    **/
    return ecc_adjlist;
}

vector<vector<NodeID>> ECC2VCC::compute_vcc_adjlist(
    const vector<vector<NodeID>> & ecc_adjlist,
    std::unordered_set<std::pair<NodeID, NodeID>> & uncovered) const {

    NodeID new_id = 1;

    std::unordered_map<std::pair<NodeID, NodeID>, NodeID> ids; 

    NodeID  last_uncovered_id = 0;

    vector<std::pair<NodeID, NodeID>> uncovered_list(uncovered.begin(), uncovered.end());
    sort(uncovered_list.begin(), uncovered_list.end()
        /**,
        [](std::pair<NodeID,NodeID> const & e1, std::pair<NodeID,NodeID> const & e2) {
            return get<0>(e1) > get<0>(e2) or (get<0>(e1) == get<0>(e2) and get<1>(e1) > get<1>(e2));
        }**/);

    std::pair<NodeID, NodeID> debug_edge;
    bool debug = false;

    for (auto edge : uncovered_list) {
        if (get<0>(edge) < get<1>(edge)) {
            last_uncovered_id = new_id;
            if (new_id == 676) {
                debug = true;
                debug_edge = edge;
                std::cout << "(" << get<0>(edge) << ", " << get<1>(edge) << ") -> " << new_id << std::endl;
            }
            ids[edge] = new_id++;
        }
    }

    // TODO: remove copy, only needed for sorting right now
    vector<vector<NodeID>> ecc_copy = ecc_adjlist;

    // TODO: start at 0... doing for compatability
    for (NodeID u = 1; u < ecc_copy.size(); u++) {
        sort(ecc_copy[u].begin(), ecc_copy[u].end());
        for (auto neighbor : ecc_copy[u]) {
            if (u < neighbor and uncovered.find(std::make_pair(u, neighbor)) == uncovered.end()) {
                if (new_id == 676) {
                    debug = true;
                    debug_edge = std::make_pair(u, neighbor);
                    std::cout << "(" << u << ", " << neighbor << ") -> " << new_id << std::endl;
                }
                ids[std::make_pair(u, neighbor)] = new_id++;
            }
        }
    }

    size_t converted_edges = 0;

    // TODO: if index by zero, remove + 1
    vector<vector<NodeID>> vcc_adjlist (last_uncovered_id + 1);

    vector<NodeID> common;
    // for each edge
    for (auto id : ids) {
        common.clear();

        // currently evaluating the edge (u,v) in input graph
        auto edge = get<0>(id);
        NodeID u = get<0>(edge);
        NodeID v = get<1>(edge);
        NodeID edge_id = get<1>(id);

        static vector<NodeID> indicator(ids.size() + 1, 0);
        for (auto neigh1 : ecc_copy[u]) indicator[neigh1] = u;
        for (auto neigh2 : ecc_copy[v]) {
            if (indicator[neigh2] == u) common.push_back(neigh2);
        }

        // Precondition: common neighborhood is sorted
        // TODO: remove this if necessary
        sort(common.begin(), common.end());

        size_t neigh_index = 0;
        for (auto neigh : common) {
            // each common neighbor indicates a triangle.
            std::pair<NodeID, NodeID> e1 = u < neigh ?
                std::make_pair(u,neigh) : std::make_pair(neigh,u);
            std::pair<NodeID, NodeID> e2 = v < neigh ?
                std::make_pair(v,neigh) : std::make_pair(neigh,v);

            assert(ids.find(e1) != ids.end() && ids.find(e2) != ids.end());
            NodeID id1 = ids[e1];
            NodeID id2 = ids[e2];

            // each of these edges are uncovered.
            if (id1 <= last_uncovered_id && id2 <= last_uncovered_id) {
                //std::cout << "Adding triangle edge (" << id1 << "," << id2 << ")" << std::endl;
                vcc_adjlist[id1-1].push_back(id2-1);
                vcc_adjlist[id2-1].push_back(id1-1);
                converted_edges += 1;
            }

            // if the current edge is covered, we're done
            if (edge_id > last_uncovered_id) continue;

            // otherwise common neighbors that are adjacent to each other 
            // form a 4-clique: must check and add.
            for (size_t i = neigh_index + 1; i < common.size(); i++) {
                std::pair<NodeID,NodeID> other_edge = std::make_pair(neigh, common[i]);
                if (ids.find(other_edge) != ids.end() && u < neigh) {
                    NodeID id1 = ids[edge];
                    NodeID id2 = ids[other_edge];
                    // if both are uncovered, then add new edge
                    if (id1 <= last_uncovered_id && id2 <= last_uncovered_id) {
                        //std::cout << "Adding 4-clique edge (" << id1 << "," << id2 << ")" << std::endl;
                        vcc_adjlist[id1-1].push_back(id2-1);
                        vcc_adjlist[id2-1].push_back(id1-1);
                        converted_edges++;
                    }
                }
            }

            neigh_index++;
        }
    }

    // for 0-based indexing
    vcc_adjlist.pop_back();

    for (vector<NodeID> & neighbors : vcc_adjlist) {
        sort(neighbors.begin(), neighbors.end());
        /**
        if (!neighbors.empty()) {
            std::cout << neighbors[0];
        }
        for (size_t i = 1; i < neighbors.size(); i++) {
            std::cout << " " << neighbors[i];
        }
        std::cout << std::endl;
        **/
    }

    std::cerr << "vcc_converted_vertices=" << vcc_adjlist.size() << std::endl;
    std::cerr << "vcc_converted_edges=" << converted_edges << std::endl;

    return vcc_adjlist;
}

vector<vector<NodeID>> ECC2VCC::ecc_to_vcc() const {

    std::unordered_set<std::pair<NodeID, NodeID>> uncovered;
    vector<vector<NodeID>> ecc_adjlist = compute_ecc_adjlist(uncovered);

    vector<vector<NodeID>> vcc_adjlist = compute_vcc_adjlist(ecc_adjlist, uncovered);

    return vcc_adjlist;
}
