// Shasta.
#include "Assembler.hpp"
#include "performanceLog.hpp"
#include "ReadLoader.hpp"
#include "timestamp.hpp"
using namespace shasta;

// Standard libraries.
#include "algorithm.hpp"
#include "iterator.hpp"


// Add reads.
// The reads are added to those already previously present.
void Assembler::addReads(
    const string& fileName,
    uint64_t minReadLength,
    bool noCache,
    const size_t threadCount)
{
    reads->checkReadsAreOpen();
    reads->checkReadNamesAreOpen();

    ReadLoader readLoader(
        fileName,
        0,  // Read representation
        minReadLength,
        noCache,
        threadCount,
        largeDataFileNamePrefix,
        largeDataPageSize,
        *reads);

    reads->checkSanity();
    reads->computeReadLengthHistogram();

    cout << "Discarded read statistics for file " << fileName << ":" << endl;
    cout << "    Discarded " << readLoader.discardedInvalidBaseReadCount <<
        " reads containing invalid bases for a total " <<
        readLoader.discardedInvalidBaseBaseCount << " valid bases." << endl;
    cout << "    Discarded " << readLoader.discardedShortReadReadCount <<
        " reads shorter than " << minReadLength <<
        " bases for a total " << readLoader.discardedShortReadBaseCount << " bases." << endl;
    cout << "    Discarded " << readLoader.discardedBadRepeatCountReadCount <<
        " reads containing repeat counts 256 or more" <<
        " for a total " << readLoader.discardedBadRepeatCountBaseCount << " bases." << endl;

    // Increment the discarded reads statistics.
    assemblerInfo->discardedInvalidBaseReadCount += readLoader.discardedInvalidBaseReadCount;
    assemblerInfo->discardedInvalidBaseBaseCount += readLoader.discardedInvalidBaseBaseCount;
    assemblerInfo->discardedShortReadReadCount += readLoader.discardedShortReadReadCount;
    assemblerInfo->discardedShortReadBaseCount += readLoader.discardedShortReadBaseCount;
    assemblerInfo->discardedBadRepeatCountReadCount += readLoader.discardedBadRepeatCountReadCount;
    assemblerInfo->discardedBadRepeatCountBaseCount += readLoader.discardedBadRepeatCountBaseCount;
    assemblerInfo->minReadLength = minReadLength;
}


// Create a histogram of read lengths.
// All lengths here are raw sequence lengths
// (length of the original read), not lengths
// in run-length representation.
void Assembler::histogramReadLength(const string& fileName)
{
    reads->computeReadLengthHistogram();
    reads->writeReadLengthHistogram(fileName);

    cout << "Discarded read statistics for all input files:" << endl;;
    cout << "    Discarded " << assemblerInfo->discardedInvalidBaseReadCount <<
        " reads containing invalid bases for a total " <<
        assemblerInfo->discardedInvalidBaseBaseCount << " valid bases." << endl;
    cout << "    Discarded " << assemblerInfo->discardedShortReadReadCount <<
        " short reads for a total " <<
        assemblerInfo->discardedShortReadBaseCount << " bases." << endl;
    cout << "    Discarded " << assemblerInfo->discardedBadRepeatCountReadCount <<
        " reads containing repeat counts 256 or more" <<
        " for a total " << assemblerInfo->discardedBadRepeatCountBaseCount << " bases." << endl;

    cout << "Read statistics for reads that will be used in this assembly:" << endl;
    cout << "    Total number of reads is " << reads->readCount() << "." << endl;
    cout << "    Total number of raw bases is " << reads->getTotalBaseCount() << "." << endl;
    cout << "    Average read length is " << double(reads->getTotalBaseCount()) / double(reads->readCount());
    cout << " bases." << endl;
    cout << "    N50 for read length is " << reads->getN50() << " bases." << endl;

    // Store read statistics in AssemblerInfo.
    assemblerInfo->readCount = reads->readCount();
    assemblerInfo->baseCount = reads->getTotalBaseCount();
    assemblerInfo->readN50 = reads->getN50();
}



uint64_t Assembler::adjustCoverageAndGetNewMinReadLength(uint64_t desiredCoverage) {
    cout << timestamp << "Adjusting for desired coverage." << endl;
    cout << "Desired Coverage: " << desiredCoverage << endl;
    uint64_t cumulativeBaseCount = reads->getTotalBaseCount();

    assemblerInfo->minReadLength = 0ULL;

    if (desiredCoverage > cumulativeBaseCount) {
        return assemblerInfo->minReadLength;
    }

    const auto& histogram = reads->getReadLengthHistogram();
    uint64_t lastLength = 0;

    for (uint64_t length = 0; length < histogram.size(); length++) {
        const uint64_t frequency = histogram[length];
        if (frequency) {
            const uint64_t baseCount = frequency * length;
            if (cumulativeBaseCount > desiredCoverage) {
                cumulativeBaseCount -= baseCount;
                lastLength = length;
                continue;
            }

            assemblerInfo->minReadLength = lastLength;
            break;
        }
    }

    cout << "Setting minReadLength to " + to_string(assemblerInfo->minReadLength) +
        " to get desired coverage." << endl;

    // Rename existing memory mapped files to avoid overwriting data.
    reads->rename();

    unique_ptr<Reads> newReads = make_unique<Reads>();
    newReads->createNew(
        0,  // Read representation
        largeDataName("Reads"),
        largeDataName("ReadNames"),
        largeDataName("ReadMetaData"),
        largeDataName("ReadFlags"),
        largeDataName("ReadIdsSortedByName"),
        largeDataPageSize
    );

    newReads->copyDataForReadsLongerThan(
        getReads(),
        assemblerInfo->minReadLength,
        assemblerInfo->discardedShortReadReadCount,
        assemblerInfo->discardedShortReadBaseCount
    );

    reads->remove();

    reads = std::move(newReads);

    // Re-compute the histogram.
    reads->computeReadLengthHistogram();

    cout << timestamp << "Done adjusting for desired coverage." << endl;

    return assemblerInfo->minReadLength;
}



void Assembler::computeReadIdsSortedByName()
{
    reads->computeReadIdsSortedByName();
}




// Find duplicate reads, as determined by name (not sequence).
// This also sets the isDuplicate and discardDueToDuplicates read flags
// and summarizes what it found Duplicates.csv.
void Assembler::findDuplicateReads(const string& handleDuplicates)
{
    reads->findDuplicates(handleDuplicates);
}
