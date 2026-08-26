#pragma once

// This finds the block structure of a rectangular matrix.
// It creates a vector of block, where each block
// is a pair of vector<uint64_t>. The two vectors
// in the pair contain the i and j indexes for the block.

// Shasta2.
#include "DisjointSets.hpp"

// Standard library.
#include "SHASTA2_ASSERT.hpp"
#include "utility.hpp"
#include "vector.hpp"

namespace shasta2 {
    template<class T> void findBlockStructure(
        const vector< vector<T> >& matrix,
        vector< pair< vector<uint64_t>, vector<uint64_t> > >& blocks
        );
}



template<class T> void shasta2::findBlockStructure(
    const vector< vector<T> >& matrix,
    vector< pair< vector<uint64_t>, vector<uint64_t> > >& blocks
    )
{
    // Get the matrix size and do sanity check.
    const uint64_t n = matrix.size();
    SHASTA2_ASSERT(n > 0);
    const uint64_t m = matrix.front().size();
    SHASTA2_ASSERT(m > 0);
    for(const auto& v: matrix) {
        SHASTA2_ASSERT(v.size() == m);
    }

    // Find connected components of a bipartite graph in which each vertex
    // represents a i or j index. A i and j index are in the
    // same connected component if the matrix has a non-zero element i,j.
    // i indexes are at positions [0,n).
    // j indexes are at positions [n, n+m).
    const T zero = T{};
    DisjointSets disjointSets(n + m);
    for(uint64_t i=0; i<n; i++) {
        for(uint64_t j=0; j<m; j++) {
            if(matrix[i][j] != zero) {
                disjointSets.unionSet(i, j + n);
            }
        }
    }
    vector< vector<uint64_t> > components;
    disjointSets.gatherComponents(1, components);

    // Each connected component generates a block.
    blocks.clear();
    for(const vector<uint64_t>& component: components) {
        auto& block = blocks.emplace_back();
        auto& iIndexes = block.first;
        auto& jIndexes = block.second;
        for(uint64_t k: component) {
            if(k < n) {
                iIndexes.push_back(k);
            } else {
                jIndexes.push_back(k - n);
            }
        }
    }

}
