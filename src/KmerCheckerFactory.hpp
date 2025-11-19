#pragma once

// Shasta.
#include "KmerChecker.hpp"
#include "memory.hpp"

namespace shasta {
    class KmerCheckerFactory;

    class KmerChecker;
    class MappedMemoryOwner;

}



// The KmerCheckerFactory knows how to create the appropriate
// type of KmerChecker for the options used.
class shasta::KmerCheckerFactory {
public:

    static shared_ptr<KmerChecker> createNew(
        uint64_t k,
        double markerDensity,
        const MappedMemoryOwner&);

    static shared_ptr<KmerChecker> createFromBinaryData(const MappedMemoryOwner&);
};


