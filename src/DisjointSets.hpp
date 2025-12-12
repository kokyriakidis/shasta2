#pragma once

// A wrapper class for boost::disjint_sets<uint64_t>.

// Boost libraries.
#include <boost/pending/disjoint_sets.hpp>

// Standard library.
#include "vector.hpp"

namespace shasta2 {
    class DisjointSets;
}



class shasta2::DisjointSets {
public:

    DisjointSets(uint64_t n);

    // This causes each set to consist of a single element.
    // It is called automatically by the constructor but
    // can be called again to reuse the same DisjointSets object.
    void initializeDisconnected();

    // These call the corresponding functions of boost::disjoint_sets.
    void unionSet(uint64_t, uint64_t);
    void link(uint64_t, uint64_t);
    uint64_t findSet(uint64_t);

    // Find all components of size at least equal to minComponentSize
    // and gather them sorted by decreasing size.
    // Each component is sorted.
    void gatherComponents(
        uint64_t minComponentSize,
        vector< vector<uint64_t> >&);

private:
    vector<uint64_t> rank;
    vector<uint64_t> parent;
    boost::disjoint_sets<uint64_t*, uint64_t*> disjointSets;
};
