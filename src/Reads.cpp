// Shasta
#include "Reads.hpp"
#include "ReadId.hpp"

// Standard Library
#include "fstream.hpp"
#include "tuple.hpp"

using namespace shasta;


void Reads::createNew(
    const string& readsDataName,
    const string& readNamesDataName,
    const string& readMetaDataDataName,
    const string& readFlagsDataName,
    const string& readIdsSortedByNameDataName,
    uint64_t largeDataPageSize)
{
    reads.createNew(readsDataName, largeDataPageSize);
    readNames.createNew(readNamesDataName, largeDataPageSize);
    readMetaData.createNew(readMetaDataDataName, largeDataPageSize);
    readFlags.createNew(readFlagsDataName, largeDataPageSize);
    readIdsSortedByName.createNew(readIdsSortedByNameDataName, largeDataPageSize);
}

void Reads::access(
    const string& readsDataName,
    const string& readNamesDataName,
    const string& readMetaDataDataName,
    const string& readFlagsDataName,
    const string& readIdsSortedByNameDataName)
{
    reads.accessExistingReadWrite(readsDataName);
    readNames.accessExistingReadWrite(readNamesDataName);
    readMetaData.accessExistingReadWrite(readMetaDataDataName);
    readFlags.accessExistingReadWrite(readFlagsDataName);
    readIdsSortedByName.accessExistingReadWrite(readIdsSortedByNameDataName);
}


void Reads::rename() {
    const string suffix = "_old";
    const string readsDataName = reads.getName();
    const string readNamesDataName = readNames.getName();
    const string readMetaDataDataName = readMetaData.getName();
    const string readFlagsDataName = readFlags.fileName;

    // No need to rename if anonymous memory mode is used.
    if (!readsDataName.empty()) {
        reads.rename(readsDataName + suffix);
    }
    if (!readNamesDataName.empty()) {
        readNames.rename(readNamesDataName + suffix);
    }
    if (!readMetaDataDataName.empty()) {
        readMetaData.rename(readMetaDataDataName + suffix);
    }
    if (!readFlagsDataName.empty()) {
        readFlags.rename(readFlagsDataName + suffix);
    }
}


void Reads::copyDataForReadsLongerThan(
    const Reads& rhs,
    uint64_t newMinReadLength,
    uint64_t& discardedShortReadCount,
    uint64_t& discardedShortReadBases
) {
    for(ReadId id = 0; id < rhs.readCount(); id++) {
        const auto len = rhs.getReadRawSequenceLength(id);
        if (len >= newMinReadLength) {
            // Copy over stuff.
            readNames.appendVector(rhs.readNames.begin(id), rhs.readNames.end(id));
            readMetaData.appendVector(rhs.readMetaData.begin(id), rhs.readMetaData.end(id));
            reads.append(rhs.reads[id]);
        } else {
            discardedShortReadCount++;
            discardedShortReadBases += len;
        }
    }

    reads.unreserve();
    readNames.unreserve();
    readMetaData.unreserve();
    readFlags.reserveAndResize(reads.size());
}



void Reads::remove() {
    reads.remove();
    readNames.remove();
    readMetaData.remove();
    readFlags.remove();
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
vector<Base> Reads::getOrientedReadRawSequence(OrientedReadId orientedReadId) const
{
    // The sequence we will return;
    vector<Base> sequence;

    // The number of bases stored.
    const uint32_t storedBaseCount = uint32_t(reads[orientedReadId.getReadId()].baseCount);

    if(representation == 0) {

        // We are storing the raw sequence of the read.
        for(uint32_t position=0; position<storedBaseCount; position++) {
            const Base base  = getOrientedReadBase(orientedReadId, position);
            sequence.push_back(base);
        }

    } else {
        SHASTA_ASSERT(0);
    }


    return sequence;
}

// Return the length of the raw sequence of a read.
// If using the run-length representation of reads, this counts each
// base a number of times equal to its repeat count.
uint64_t Reads::getReadRawSequenceLength(ReadId readId) const
{
    return reads[readId].baseCount;
}



// Function to write one or all reads in Fasta format.
void Reads::writeReads(const string& fileName)
{
    ofstream file(fileName);
    for(ReadId readId=0; readId<reads.size(); readId++) {
        writeRead(readId, file);
    }
}


void Reads::writeRead(ReadId readId, const string& fileName)
{
    ofstream file(fileName);
    writeRead(readId, file);
}


void Reads::writeRead(ReadId readId, ostream& file)
{
    checkReadsAreOpen();
    checkReadNamesAreOpen();
    checkReadId(readId);

    const vector<Base> rawSequence = getOrientedReadRawSequence(OrientedReadId(readId, 0));
    const auto readName = readNames[readId];
    const auto metaData = readMetaData[readId];

    file << ">";
    copy(readName.begin(), readName.end(), ostream_iterator<char>(file));
    file << " " << readId;
    file << " " << rawSequence.size();
    if(metaData.size() > 0) {
        file << " ";
        copy(metaData.begin(), metaData.end(), ostream_iterator<char>(file));
    }
    file << "\n";
    copy(rawSequence.begin(), rawSequence.end(), ostream_iterator<Base>(file));
    file << "\n";

}


void Reads::writeOrientedRead(ReadId readId, Strand strand, const string& fileName)
{
    writeOrientedRead(OrientedReadId(readId, strand), fileName);
}


void Reads::writeOrientedRead(OrientedReadId orientedReadId, const string& fileName)
{
    ofstream file(fileName);
    writeOrientedRead(orientedReadId, file);
}


void Reads::writeOrientedRead(OrientedReadId orientedReadId, ostream& file)
{
    checkReadsAreOpen();
    checkReadNamesAreOpen();

    const vector<Base> rawSequence = getOrientedReadRawSequence(orientedReadId);
    const auto readName = readNames[orientedReadId.getReadId()];

    file << ">" << orientedReadId;
    file << " " << rawSequence.size() << " ";
    copy(readName.begin(), readName.end(), ostream_iterator<char>(file));
    file << "\n";
    copy(rawSequence.begin(), rawSequence.end(), ostream_iterator<Base>(file));
    file << "\n";
}

// Create a histogram of read lengths.
// All lengths here are raw sequence lengths
// (length of the original read), not lengths
// in run-length representation.
void Reads::computeReadLengthHistogram() {
    checkReadsAreOpen();
    histogram.clear();

    // Create the histogram.
    // It contains the number of reads of each length.
    // Indexed by the length.
    const ReadId totalReadCount = readCount();
    totalBaseCount = 0;
    
    for(ReadId readId=0; readId<totalReadCount; readId++) {
        const uint64_t length = getReadRawSequenceLength(readId);
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

