#pragma once

#include "Detangler.hpp"

#include "cstdint.hpp"

namespace shasta {
    class PermutationDetangler;
}


class shasta::PermutationDetangler : public Detangler {
public:
    bool operator()(Tangle&);

    PermutationDetangler(uint64_t minCommonCoverage) :
        minCommonCoverage(minCommonCoverage) {}
    uint64_t minCommonCoverage;
};
