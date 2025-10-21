// Shasta
#include "Reads.hpp"
#include "extractKmer128.hpp"
#include "ReadId.hpp"

// Standard Library
#include "fstream.hpp"
#include "tuple.hpp"

using namespace shasta;


void Reads::createNew(
    const string& readsDataName,
    const string& readNamesDataName,
    const string& readIdsSortedByNameDataName,
    uint64_t largeDataPageSize)
{
    reads.createNew(readsDataName, largeDataPageSize);
    readNames.createNew(readNamesDataName, largeDataPageSize);
    readIdsSortedByName.createNew(readIdsSortedByNameDataName, largeDataPageSize);
}

void Reads::access(
    const string& readsDataName,
    const string& readNamesDataName,
    const string& readIdsSortedByNameDataName)
{
    reads.accessExistingReadOnly(readsDataName);
    readNames.accessExistingReadOnly(readNamesDataName);
    readIdsSortedByName.accessExistingReadOnly(readIdsSortedByNameDataName);
}



void Reads::remove() {
    reads.remove();
    readNames.remove();
    readIdsSortedByName.remove();
}



// Return a base of an oriented read.
Base Reads::getOrientedReadBase(
    OrientedReadId orientedReadId,
    uint32_t position) const
{
    const auto& read = reads[orientedReadId.getReadId()];
    if(orientedReadId.getStrand() == 0) {
        return read[position];
    } else {
        return read[read.baseCount-1-position].complement();
    }
}



// Return a vector containing the raw sequence of an oriented read.
vector<Base> Reads::getOrientedReadSequence(OrientedReadId orientedReadId) const
{
    // The sequence we will return;
    vector<Base> sequence;

    const uint32_t baseCount = uint32_t(reads[orientedReadId.getReadId()].baseCount);
    sequence.reserve(baseCount);

    for(uint32_t position=0; position < baseCount; position++) {
        const Base base  = getOrientedReadBase(orientedReadId, position);
        sequence.push_back(base);
    }

    return sequence;
}



uint64_t Reads::getReadSequenceLength(ReadId readId) const
{
    return reads[readId].baseCount;
}



// Create a histogram of read lengths.
void Reads::computeReadLengthHistogram() {
    checkReadsAreOpen();
    histogram.clear();

    // Create the histogram.
    // It contains the number of reads of each length.
    // Indexed by the length.
    const ReadId totalReadCount = readCount();
    totalBaseCount = 0;
    
    for(ReadId readId=0; readId<totalReadCount; readId++) {
        const uint64_t length = getReadSequenceLength(readId);
        totalBaseCount += length;
        if(histogram.size() <= length) {
            histogram.resize(length+1, 0);
        }
        ++(histogram[length]);
    }

    // Binned histogram
    binnedHistogram.clear();
    const uint64_t binWidth = 1000;
    for(uint64_t length=0; length<histogram.size(); length++) {
        const uint64_t readCount = histogram[length];
        if(readCount) {
            const uint64_t bin = length / binWidth;
            if(binnedHistogram.size() <= bin) {
                binnedHistogram.resize(bin+1, make_pair(0, 0));
            }
            binnedHistogram[bin].first += readCount;
            binnedHistogram[bin].second += readCount * length;
        }
    }
}

void Reads::writeReadLengthHistogram(const string& fileName) {
    checkReadsAreOpen();
    const ReadId totalReadCount = readCount();
    
    n50 = 0;
    {
        ofstream csv(fileName);
        csv << "Length,Reads,Bases,CumulativeReads,CumulativeBases,"
            "FractionalCumulativeReads,FractionalCumulativeBases,\n";
        uint64_t cumulativeReadCount = totalReadCount;
        uint64_t cumulativeBaseCount = totalBaseCount;
        for(uint64_t length=0; length<histogram.size(); length++) {
            const uint64_t frequency = histogram[length];
            if(frequency) {
                const  uint64_t baseCount = frequency * length;
                const double cumulativeReadFraction =
                    double(cumulativeReadCount)/double(totalReadCount);
                const double comulativeBaseFraction =
                    double(cumulativeBaseCount)/double(totalBaseCount);
                csv << length << "," << frequency << "," << baseCount << ",";
                csv << cumulativeReadCount << "," << cumulativeBaseCount << ",";
                csv << cumulativeReadFraction << ",";
                csv << comulativeBaseFraction << "\n";
                cumulativeReadCount -= frequency;
                cumulativeBaseCount -= baseCount;
                if(comulativeBaseFraction > 0.5) {
                    n50 = length;
                }
            }
        }
        
        SHASTA_ASSERT(cumulativeReadCount == 0);
        SHASTA_ASSERT(cumulativeBaseCount == 0);
    }

    // Binned Histogram.
    {
        const int binWidth = 1000;
        ofstream csv("Binned-" + fileName);
        csv << "LengthBegin,LengthEnd,Reads,Bases,CumulativeReads,CumulativeBases,"
            "FractionalCumulativeReads,FractionalCumulativeBases,\n";
        uint64_t cumulativeReadCount = totalReadCount;
        uint64_t cumulativeBaseCount = totalBaseCount;
        for(uint64_t bin=0; bin<binnedHistogram.size(); bin++) {
            const auto& histogramBin = binnedHistogram[bin];
            const uint64_t readCount = histogramBin.first;
            const uint64_t baseCount = histogramBin.second;
            const double cumulativeReadFraction =
                double(cumulativeReadCount)/double(totalReadCount);
            const double comulativeBaseFraction =
                double(cumulativeBaseCount)/double(totalBaseCount);
            csv << bin*binWidth << ",";
            csv << (bin+1)*binWidth << ",";
            csv << readCount << "," << baseCount << ",";
            csv << cumulativeReadCount << "," << cumulativeBaseCount << ",";
            csv << cumulativeReadFraction << ",";
            csv << comulativeBaseFraction << "\n";
            cumulativeReadCount -= readCount;
            cumulativeBaseCount -= baseCount;
        }
        SHASTA_ASSERT(cumulativeReadCount == 0);
        SHASTA_ASSERT(cumulativeBaseCount == 0);
    }

}



void Reads::computeReadIdsSortedByName()
{
    // Store ReadId's in numerical order.
    readIdsSortedByName.resize(readCount());
    for(ReadId readId=0; readId<readCount(); readId++) {
        readIdsSortedByName[readId] = readId;
    }

    // Sort them by name.
    sort(readIdsSortedByName.begin(), readIdsSortedByName.end(),
        OrderReadsByName(readNames));

    // Check for duplicate names.
    for(ReadId readId=1; readId<readCount(); readId++) {
        if(getReadName(readId - 1) == getReadName(readId)) {
            string message = "Duplicate read name: ";
            const auto name = getReadName(readId);
            copy(name.begin(), name.end(), back_inserter(message));
            throw runtime_error(message);
        }
    }
}



// Get a ReadId given a read name.
// This uses a binary search in readIdsSortedByName.
ReadId Reads::getReadId(const string& readName) const
{
    const auto begin = readName.data();
    const auto end = begin + readName.size();
    span<const char> s(begin, end);
    return getReadId(s);
}
ReadId Reads::getReadId(const span<const char>& readName) const
{
    const auto begin = readIdsSortedByName.begin();
    const auto end = readIdsSortedByName.end();
    auto it = std::lower_bound(begin, end, readName, OrderReadsByName(readNames));
    if(it == end) {
        return invalidReadId;
    }
    const ReadId readId = *it;
    if(readNames[readId] == readName) {
        return readId;
    } else {
        return invalidReadId;
    }
}



Kmer Reads::getKmer(uint64_t k, OrientedReadId orientedReadId, uint32_t position) const
{
    // Get the ReadId and check that it is valid.
    const ReadId readId = orientedReadId.getReadId();
    if(readId >= readCount()) {
        throw runtime_error("Invalid ReadId.");
    }

    // Get the read sequence and check that the position is valid.
    const LongBaseSequenceView read = getRead(readId);
    const uint64_t readLength = read.baseCount;
    if(position + k > readLength) {
        throw runtime_error("Invalid position in read.");
    }

    // Extract the Kmer.
    const uint64_t strand = orientedReadId.getStrand();
    if(strand == 0) {
        Kmer kmer;
        extractKmer128(read, position, k, kmer);
        return kmer;
    } else {
        Kmer kmer;
        extractKmer128(read, readLength - 1 - position, k, kmer);
        kmer.reverseComplement(k);
        return kmer;
    }
}
