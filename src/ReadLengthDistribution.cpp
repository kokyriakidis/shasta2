// Shasta.
#include "ReadLengthDistribution.hpp"
#include "Reads.hpp"
#include "shastaTypes.hpp"
using namespace shasta;

// Standard library.
#include "fstream.hpp"


// Initial creation.
ReadLengthDistribution::ReadLengthDistribution(
    const Reads& reads,
    const MappedMemoryOwner& mappedMemoryOwner):
    MappedMemoryOwner(mappedMemoryOwner)
{
    vector<uint64_t> histogram;
    computeReadLengthHistogram(reads, histogram);
    const uint64_t binCount = histogram. size() + 1;

    data.createNew(largeDataName("ReadLengthDistribution"), largeDataPageSize);
    data.resize(binCount);

    // Fill in the P values.
    const uint64_t readCount = reads.readCount();
    for(uint64_t bin=0; bin<binCount; bin++) {
        data[bin].P = double(histogram[bin]) / double(readCount);
    }

    // Fill in the Q and R values.
    data.back().Q = 0.;
    data.back().R = 0.;
    for(uint64_t bin=binCount-1; /* Check later */; bin--) {
        const uint64_t lengthBegin = bin * binWidth;
        const uint64_t lengthEnd = (bin + 1) * binWidth;
        const double averageLength = double(lengthBegin + lengthEnd) / 2.;

        data[bin].Q = data[bin+1].Q + data[bin].P;
        data[bin].R = data[bin+1].R + data[bin].P / averageLength;

        // Must check here because loop variable is unsigned.
        if(bin == 0) {
            break;
        }
    }

    // Compute the coverage correlation values.
    for(uint64_t bin=0; bin < binCount; bin++) {
        const uint64_t lengthBegin = bin * binWidth;
        auto& entry = data[bin];
        entry.coverageCorrelation = entry.Q - double(lengthBegin) * entry.R;
    }


    write(histogram);
}





// Access existing from binary data.
ReadLengthDistribution::ReadLengthDistribution(
    const MappedMemoryOwner& mappedMemoryOwner) :
    MappedMemoryOwner(mappedMemoryOwner)
{
    data.accessExistingReadOnly(largeDataName("ReadLengthDistribution"));
}



void ReadLengthDistribution::computeReadLengthHistogram(
    const Reads& reads,
    vector<uint64_t>& histogram) const
{
    histogram.clear();

    const uint64_t readCount = reads.readCount();
    for(ReadId readId=0; readId<readCount; readId++) {
        const uint64_t length = reads.getReadSequenceLength(readId);
        const uint64_t bin = length / binWidth;

        if(histogram.size() <= bin) {
            histogram.resize(bin + 1, 0);
        }
        ++histogram[bin];
    }

}



void ReadLengthDistribution::write(const vector<uint64_t>& histogram) const
{
    ofstream csv("ReadLengthDistribution.csv");
    csv << "LengthBegin,LengthEnd,Reads,P,Q,R,Coverage correlation\n";

    for(uint64_t bin=0; bin<data.size(); bin++) {
        const uint64_t frequency = histogram[bin];

        const uint64_t lengthBegin = bin * binWidth;
        const uint64_t lengthEnd = (bin + 1) * binWidth;
        csv << lengthBegin << ",";
        csv << lengthEnd << ",";
        csv << frequency << ",";

        csv << data[bin].P << ",";
        csv << data[bin].Q << ",";
        csv << data[bin].R << ",";
        csv << data[bin].coverageCorrelation << ",\n";
    }

}
