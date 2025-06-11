#pragma once

#include "Detangler.hpp"

#include "cstdint.hpp"

namespace shasta {
    class LikelihoodRatioDetangler;
}



class shasta::LikelihoodRatioDetangler : public Detangler {
public:
    bool operator()(Tangle&, bool doDetangle);

    LikelihoodRatioDetangler(
        uint64_t minCommonCoverage,
        const double epsilon,
        const double maxLogP,
        const double minLogPDelta);

    uint64_t minCommonCoverage;
    double epsilon;
    double maxLogP;
    double minLogPDelta;
};
