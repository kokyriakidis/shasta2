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

    SimpleDetangler(
        uint64_t minCommonCoverage,
        uint64_t minDetangleCoverage) :
        minCommonCoverage(minCommonCoverage),
        minDetangleCoverage(minDetangleCoverage)
    {}

    uint64_t minCommonCoverage;
    uint64_t minDetangleCoverage;
};
