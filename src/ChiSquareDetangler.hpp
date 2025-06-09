#pragma once

#include "Detangler.hpp"

#include "cstdint.hpp"

namespace shasta {
    class ChiSquareDetangler;
}



class shasta::ChiSquareDetangler : public Detangler {
public:
    bool operator()(Tangle3&);

    ChiSquareDetangler(
        uint64_t minCommonCoverage,
        const double epsilon,
        const double maxLogP,
        const double minLogPDelta);

    uint64_t minCommonCoverage;
    double epsilon;
    double maxLogP;
    double minLogPDelta;
};
