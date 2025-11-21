#pragma once

#include "MappedMemoryOwner.hpp"
#include "MemoryMappedVector.hpp"

namespace shasta2 {
    class ReadLengthDistribution;
    class ReadLengthDistributionEntry;

    class Reads;
}


// There is one ReadLengthDistributionEntry for each length bin.
class shasta2::ReadLengthDistributionEntry {
public:

    // P(l), the fraction of reads with length in this length bin.
    // This is a fraction of number of reads, not number of bases.
    // So it can diverge at short l.
    double P;

    // The cumulative distribution (integral of P from here to infinity).
    double Q;

    // The correlation factor (integral of P/l from here to infinity).
    double R;

    // The coverage correlation.
    // This is the probability that if a read covers position x
    // it also covers position x + delta.
    double coverageCorrelation;

};



class shasta2::ReadLengthDistribution : public MappedMemoryOwner {
public:

    // The width in bases of each length bin.
    const uint64_t binWidth = 1000;

    // A vector of ReadLengthDistributionEntry, one for each bin.
    MemoryMapped::Vector<ReadLengthDistributionEntry> data;

    // Initial creation.
    ReadLengthDistribution(const Reads&, const MappedMemoryOwner&);

    // Access existing from binary data.
    ReadLengthDistribution(const MappedMemoryOwner&);

public:
    void computeReadLengthHistogram(const Reads&, vector<uint64_t>& histogram) const;
    void write(const vector<uint64_t>& histogram) const;
};

