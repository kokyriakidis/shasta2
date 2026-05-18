#pragma once

// A KmerChecker that accepts all k-mers specified externally as markers.
// Use this when you want to supply your own set of anchor k-mers
// instead of using the hash-based random marker selection.

#include "KmerChecker.hpp"
#include "MurmurHash2.hpp"

#include <vector>

namespace shasta2 {
    class DinaraKmerChecker;
}



class shasta2::DinaraKmerChecker :
    public KmerChecker {
public:
    bool isMarker(const Kmer&) const;

    // Construct from a vector of k-mers.
    // Both each k-mer and its reverse complement are added to the table,
    // so the KmerChecker invariant (marker iff reverse complement is also a marker)
    // is guaranteed.
    DinaraKmerChecker(uint64_t k, const std::vector<Kmer>& kmers);

private:
    uint64_t k;

    // Hash table of marker k-mers, using the same MurmurHash + power-of-2
    // bucket masking pattern used elsewhere in shasta2 (see MarkerKmers).
    std::vector< std::vector<Kmer> > buckets;
    uint64_t mask;

    uint64_t findBucket(const Kmer& kmer) const
    {
        const uint64_t hashValue = MurmurHash64A(&kmer, sizeof(kmer), 1241);
        return hashValue & mask;
    }

    bool find(const Kmer&) const;
};
