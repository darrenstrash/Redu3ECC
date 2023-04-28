#pragma once
#include <vector>
#include <unordered_set>

// for hashing
#include <boost/functional/hash.hpp>

#include "graph_io.h"
#include "cover.hpp"
#include "graph.hpp"

class ECC2VCC {
    private:
        const ECCGraph & graph;
        const Cover    & cover;
    public:
        ECC2VCC(const ECCGraph & graph, const Cover & cover)
        : graph(graph)
        , cover(cover)
        { }

    std::vector<std::vector<NodeID>> ecc_to_vcc() const;

    std::vector<std::vector<NodeID>> compute_ecc_adjlist(std::unordered_set<std::pair<NodeID, NodeID>> & uncovered) const;

    std::vector<std::vector<NodeID>> compute_vcc_adjlist(std::vector<std::vector<NodeID>> const & ecc_adjlist, std::unordered_set<std::pair<NodeID, NodeID>> & uncovered) const;
};
