#pragma once

#include "Detangler.hpp"

#include "cstdint.hpp"

namespace shasta2 {
    class LikelihoodRatioDetangler;
}



class shasta2::LikelihoodRatioDetangler : public Detangler {
public:
    bool operator()(Tangle1&);

    LikelihoodRatioDetangler(
        const double epsilon,
        const double maxLogP,
        const double minLogPDelta,
        uint64_t detangleMinCoverage,
        bool onlyConsiderInjective,
        bool onlyConsiderPermutation,
        bool requireInjective,
        bool requirePermutation);

    double epsilon;
    double maxLogP;
    double minLogPDelta;
    uint64_t detangleMinCoverage;

    // If this is set, the GTest only considers hypotheses that are
    // forward injective and/or backward injective.
    bool onlyConsiderInjective;

    // If this is set, the GTest only considers hypotheses that are
    // forward injective and backward injective - that is, they are
    // permutations.
    bool onlyConsiderPermutation;

    // If this is set, all hypotheses are computed, but detangling
    // is not done if the top hypothesis is not forward injective
    // or backward injective.
    bool requireInjective = true;

    // If this is set, all hypotheses are computed, but detangling
    // is not done if the top hypothesis is not forward injective
    // and backward injective (that is, is a permutation).
    bool requirePermutation = false;
};
