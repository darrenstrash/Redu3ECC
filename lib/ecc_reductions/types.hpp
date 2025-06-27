#ifndef ECC_REDUCTIONS_TYPES
#define ECC_REDUCTIONS_TYPES

#include <unordered_map>
#include <unordered_set>
#include <cstdint>
#include <vector>
#include <algorithm>

// for hashing
#include <boost/functional/hash.hpp>

// graph types
typedef uint32_t node_t; // I might eventually need more than 4 billion nodes...

// types to clean things up

// allow hashing with node pairs
typedef std::pair<node_t,node_t> edge_t;

namespace std
{
    template<> struct hash<std::pair<node_t,node_t>>
    {
        std::size_t operator()(std::pair<node_t,node_t> const & node_pair) const noexcept
        {
            std::size_t seed = 0;
            boost::hash_combine(seed, std::get<0>(node_pair));
            boost::hash_combine(seed, std::get<1>(node_pair));
            return seed;
        }
    };
};

typedef std::vector<std::vector<node_t>> adjlist_t;

// set types
typedef std::unordered_set<edge_t> edgeset_t;
typedef std::unordered_set<node_t> nodeset_t;

// map types
typedef std::unordered_map<node_t,node_t> nodemap_t;
typedef std::unordered_map<node_t,edge_t> nodeedgemap_t;
typedef std::unordered_map<edge_t,node_t> edgenodemap_t;

typedef std::vector<node_t> clique_t;
typedef std::vector<clique_t> cover_t;

class NodeVector : public std::vector<node_t> {
public:
    // please don't call this
    // TODO: Remove calls for speed, but for now it's only slow not incorrect.
    bool contains(node_t const& node) const {
        //assert(0);
        return std::find(cbegin(), cend(), node) != cend();
    }

    void insert(node_t const& node) {
        push_back(node);
    }

    // please don't call this either
    void erase(node_t const& node) {
        assert(0);
        // two full passes. Should just perform one pass, swap to end and remove.
        auto it = std::find(cbegin(), cend(), node);
        if (it != cend()) {
            std::vector<node_t>& vec = *this;
            vec.erase(it);
        }
    }
};

typedef NodeVector node_container_t;
typedef std::unordered_map<node_t, node_container_t> adj_list_data_t;



#endif //ECC_REDUCTIONS_TYPES
