#pragma once

#include "KmerChecker.hpp"

namespace shasta2 {
    class HashedKmerChecker;
}


// The new implementation of the KmerChecker is not table based
// and uses hashing instead.
// It only supports marker generation method 0 (random generation)
// but allow marker lengths k<32.
class shasta2::HashedKmerChecker :
    public KmerChecker {
public:
    bool isMarker(const Kmer&) const;

    // Initial creation.
    HashedKmerChecker(uint64_t k, double markerDensity);

private:
    uint64_t k;
    uint32_t hashThreshold;

};


