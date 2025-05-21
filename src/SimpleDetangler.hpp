#pragma once

#include "Detangler.hpp"

#include "cstdint.hpp"

namespace shasta {
    class SimpleDetangler;
}


class shasta::SimpleDetangler : public Detangler {
public:
    bool operator()(Tangle&);
    bool operator()(Tangle2&);
    bool operator()(Tangle3&);

    // Tangle matrix elements <= detangleLowThreshold are considered insignificant.
    // Tangle matrix elements >0 detangleHighThreshold are considered significant.
    // Tangle matrix elements in-between are considered ambiguous and cause the
    // detangling to fail.
    // If there are no ambiguous elements, detangling is successful and a connection
    // is generated for each significant tangle matrix element.
    SimpleDetangler(
        uint64_t minCommonCoverage,
        uint64_t detangleLowThreshold,
        uint64_t detangleHighThreshold) :
        minCommonCoverage(minCommonCoverage),
        detangleLowThreshold(detangleLowThreshold),
        detangleHighThreshold(detangleHighThreshold)
    {}

    uint64_t minCommonCoverage;
    uint64_t detangleLowThreshold;
    uint64_t detangleHighThreshold;
};
