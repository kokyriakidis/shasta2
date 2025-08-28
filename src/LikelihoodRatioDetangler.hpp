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
        const double minLogPDelta);

    double epsilon;
    double maxLogP;
    double minLogPDelta;
};
