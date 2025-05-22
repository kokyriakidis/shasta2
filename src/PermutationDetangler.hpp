#pragma once

#include "Detangler.hpp"

#include "cstdint.hpp"

namespace shasta {
    class PermutationDetangler;
}


class shasta::PermutationDetangler : public Detangler {
public:
    bool operator()(Tangle&);
    bool operator()(Tangle3&);

    PermutationDetangler(
        uint64_t minCommonCoverage,
        const double epsilon,
        const double maxLogP,
        const double minLogPDelta);

    uint64_t minCommonCoverage;
    double epsilon;
    double maxLogP;
    double minLogPDelta;
};
