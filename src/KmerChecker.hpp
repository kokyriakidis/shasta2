#ifndef SHASTA_KMER_CHECKER_HPP
#define SHASTA_KMER_CHECKER_HPP

// Shasta.
#include "Kmer.hpp"

namespace shasta {
    class KmerChecker;
}



// The KmerChecker is an abstract class that knows how to find
// out if a k-mer is a marker.
// All implementations must guarantee that if a Kmer is a marker
// its reverse complement is also a marker.
class shasta::KmerChecker {
public:
    virtual bool isMarker(const Kmer&) const = 0;
};

#endif
