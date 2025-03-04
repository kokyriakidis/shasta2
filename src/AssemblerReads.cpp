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
    const size_t threadCount)
{
    reads->checkReadsAreOpen();
    reads->checkReadNamesAreOpen();

    ReadLoader readLoader(
        fileName,
        minReadLength,
        threadCount,
        largeDataFileNamePrefix,
        largeDataPageSize,
        *reads);

    reads->checkSanity();
    reads->computeReadLengthHistogram();
}


// Create a histogram of read lengths.
// All lengths here are raw sequence lengths
// (length of the original read), not lengths
// in run-length representation.
void Assembler::histogramReadLength(const string& fileName)
{
    reads->computeReadLengthHistogram();
    reads->writeReadLengthHistogram(fileName);


    cout << "Read statistics for reads that will be used in this assembly:" << endl;
    cout << "    Total number of reads is " << reads->readCount() << "." << endl;
    cout << "    Total number of bases is " << reads->getTotalBaseCount() << "." << endl;
    cout << "    Average read length is " << double(reads->getTotalBaseCount()) / double(reads->readCount());
    cout << " bases." << endl;
    cout << "    N50 for read length is " << reads->getN50() << " bases." << endl;

    // Store read statistics in AssemblerInfo.
    assemblerInfo->readCount = reads->readCount();
    assemblerInfo->baseCount = reads->getTotalBaseCount();
    assemblerInfo->readN50 = reads->getN50();
}



void Assembler::computeReadIdsSortedByName()
{
    reads->computeReadIdsSortedByName();
}
