// Shasta.
#include "DinaraKmerChecker.hpp"
using namespace shasta2;

#include <algorithm>
#include <bit>



DinaraKmerChecker::DinaraKmerChecker(
    uint64_t k,
    const std::vector<Kmer>& kmers) :
    k(k)
{
    // Choose the number of buckets as the smallest power of 2
    // that is >= the number of k-mers (including reverse complements).
    // This keeps the average bucket size around 1.
    const uint64_t n = 2 * kmers.size();
    const uint64_t bucketCount = std::max(uint64_t(1), std::bit_ceil(n));
    mask = bucketCount - 1;
    buckets.resize(bucketCount);

    // Insert each k-mer and its reverse complement.
    for(const Kmer& kmer: kmers) {
        buckets[findBucket(kmer)].push_back(kmer);
        const Kmer kmerRc = kmer.reverseComplement(k);
        buckets[findBucket(kmerRc)].push_back(kmerRc);
    }
}



bool DinaraKmerChecker::find(const Kmer& kmer) const
{
    const std::vector<Kmer>& bucket = buckets[findBucket(kmer)];
    return std::find(bucket.begin(), bucket.end(), kmer) != bucket.end();
}



// A k-mer is a marker if it or its reverse complement is in the table.
bool DinaraKmerChecker::isMarker(const Kmer& kmer) const
{
    if(find(kmer)) {
        return true;
    }
    return find(kmer.reverseComplement(k));
}
