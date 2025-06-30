#pragma once

#include "Detangler.hpp"
#include "SHASTA_ASSERT.hpp"

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
    // If there there are no ambiguous elements and at least one insignificant element,
    // detangling is successful and a connection is generated for each significant tangle matrix element.
    SimpleDetangler(
        uint64_t detangleLowThreshold,
        uint64_t detangleHighThreshold) :
        detangleLowThreshold(detangleLowThreshold),
        detangleHighThreshold(detangleHighThreshold)
    {
        SHASTA_ASSERT(detangleHighThreshold > detangleLowThreshold);
    }

    uint64_t detangleLowThreshold;
    uint64_t detangleHighThreshold;
};
