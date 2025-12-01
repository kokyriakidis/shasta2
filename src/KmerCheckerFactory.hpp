#pragma once

// Shasta.
#include "KmerChecker.hpp"
#include "memory.hpp"

namespace shasta2 {
    class KmerCheckerFactory;

    class KmerChecker;
}



// The KmerCheckerFactory knows how to create the appropriate
// type of KmerChecker for the options used.
class shasta2::KmerCheckerFactory {
public:

    static shared_ptr<KmerChecker> createNew(
        uint64_t k,
        double markerDensity);

};


