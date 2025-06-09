#pragma once

#include "Detangler.hpp"

#include "cstdint.hpp"

namespace shasta {
    class SimpleDetangler;
}


class shasta::SimpleDetangler : public Detangler {
public:
    bool operator()(Tangle&);

    // Tangle matrix elements <= detangleLowThreshold are considered insignificant.
    // Tangle matrix elements >= detangleHighThreshold are considered significant.
    // Tangle matrix elements in-between are considered ambiguous and cause the
    // detangling to fail.
    // If there are no ambiguous elements, detangling is successful and a connection
    // is generated for each significant tangle matrix element.
    // A maximum base offset is also enforced. If any of the significant bridge anchor pairs
    // have an offset longer than this, the detangling is not done.
    SimpleDetangler(
        uint64_t minCommonCoverage,
        uint64_t detangleLowThreshold,
        uint64_t detangleHighThreshold,
        uint64_t maxBaseOffset) :
        minCommonCoverage(minCommonCoverage),
        detangleLowThreshold(detangleLowThreshold),
        detangleHighThreshold(detangleHighThreshold),
        maxBaseOffset(maxBaseOffset)
    {}

    uint64_t minCommonCoverage;
    uint64_t detangleLowThreshold;
    uint64_t detangleHighThreshold;
    uint64_t maxBaseOffset;
};
