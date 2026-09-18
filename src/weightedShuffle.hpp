#pragma once

// Shasta2
#include "orderPairs.hpp"

// Standard library.
#include "algorithm.hpp"
#include "cstdint.hpp"
#include <random>
#include "utility.hpp"
#include "vector.hpp"

namespace shasta2 {
    template<class Generator> void weightedShuffle(
        const vector<double>& weights,
        Generator&,
        vector<uint64_t>& shuffle);

    void testWeightedShuffle();
}



// Weighted shuffle using the algorithm by Efraimidis and Spirakis, 2005
// https://utopia.duth.gr/~pefraimi/research/data/2007EncOfAlg.pdf
// Input:
// - The weights of the items to be shuffled.
// - A random bit generator (for example std::mt19937).
// Output:
// - The indexes of the elements in the weight vectors,
//   as prescribed by the shuffle.

template<class Generator> void shasta2::weightedShuffle(
    const vector<double>& weights,
    Generator& generator,
    vector<uint64_t>& shuffle)
{
    std::uniform_real_distribution<double> distribution;
    const uint64_t n = weights.size();

    // Gather pairs (index, key) where key=random^(1/weight).
    vector< pair<uint64_t, double> > indexesWithKeys(n);
    for(uint64_t i=0; i<n; i++) {
        const double weight = weights[i];
        const double u = distribution(generator);
        const double key = std::pow(u, 1. / weight);
        indexesWithKeys[i] = {i, key};
    }

    // Sort by key.
    sort(indexesWithKeys.begin(), indexesWithKeys.end(),
        OrderPairsBySecondOnlyGreater<uint64_t, double>());

    // Copy the indexes in this order.
    shuffle.resize(n);
    for(uint64_t i=0; i<n; i++) {
        shuffle[i] = indexesWithKeys[i].first;
    }

}

