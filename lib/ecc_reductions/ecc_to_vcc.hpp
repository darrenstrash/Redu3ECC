#pragma once
#include <vector>
#include <unordered_set>
#include <unordered_map>

#include "graph_io.h"
#include "cover.hpp"
#include "graph.hpp"

typedef std::vector<NodeID> clique_t;

class ECC2VCC {
    private:
        const ECCGraph & graph;
        const Cover    & cover;
    public:
        ECC2VCC(const ECCGraph & graph, const Cover & cover)
        : graph(graph)
        , cover(cover)
        { }

    // return vcc adjlist if we don't need original cover
    adjlist_t ecc_to_vcc() const;

    adjlist_t ecc_to_vcc(nodeedgemap_t & vcc_vertex_to_ecc_edge_map) const;

    adjlist_t compute_ecc_adjlist(
            edgeset_t & uncovered,
            nodemap_t & to_old_id) const;

    adjlist_t compute_vcc_adjlist(
            adjlist_t const & ecc_adjlist, 
            edgeset_t const & uncovered,
            nodeedgemap_t   & vertex_to_edge_map) const;

    void add_vcc_cliques_to_ecc_cover(
            std::vector<clique_t> const & vcc_cliques,
            nodeedgemap_t & v_to_e_map);
};
