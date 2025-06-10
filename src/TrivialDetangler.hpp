#pragma once

#include "Detangler.hpp"

#include "cstdint.hpp"

namespace shasta {
    class TrivialDetangler;
}


class shasta::TrivialDetangler : public Detangler {
public:
    bool operator()(Tangle&, bool doDetangle);

    TrivialDetangler(uint64_t minCommonCoverage) :
        minCommonCoverage(minCommonCoverage) {}
    uint64_t minCommonCoverage;
};
