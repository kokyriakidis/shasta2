#pragma once

#include "Detangler.hpp"

#include "cstdint.hpp"

namespace shasta {
    class LikelihoodRatioDetangler;
}



class shasta::LikelihoodRatioDetangler : public Detangler {
public:
    bool operator()(Tangle1&);

    LikelihoodRatioDetangler(
        const double epsilon,
        const double maxLogP,
        const double minLogPDelta,
        uint64_t detangleMinCoverage,
        bool requireInjective,
        bool requirePermutation);

    double epsilon;
    double maxLogP;
    double minLogPDelta;
    uint64_t detangleMinCoverage;

    // If this is set, all hypotheses are computed, but detangling
    // is not done if the top hypothesis is not forward injective
    // or backward injective.
    bool requireInjective = true;

    // If this is set, all hypotheses are computed, but detangling
    // is not done if the top hypothesis is not forward injective
    // and backward injective (that is, is a permutation).
    bool requirePermutation = false;
};
