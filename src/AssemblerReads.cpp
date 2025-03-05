// Shasta.
#include "Assembler.hpp"
#include "performanceLog.hpp"
#include "ReadLoader.hpp"
#include "timestamp.hpp"
using namespace shasta;

// Standard libraries.
#include "algorithm.hpp"
#include <cmath>
#include "iterator.hpp"


// Add reads.
// The reads are added to those already previously present.
void Assembler::addReads(
    const string& fileName,
    uint64_t minReadLength,
    const size_t threadCount)
{
    readsPointer->checkReadsAreOpen();
    readsPointer->checkReadNamesAreOpen();

    ReadLoader readLoader(
        fileName,
        minReadLength,
        threadCount,
        largeDataFileNamePrefix,
        largeDataPageSize,
        *readsPointer);

    readsPointer->checkSanity();
    readsPointer->computeReadLengthHistogram();
}


// Create a histogram of read lengths.
// All lengths here are raw sequence lengths
// (length of the original read), not lengths
// in run-length representation.
void Assembler::histogramReadLength(const string& fileName)
{
    readsPointer->computeReadLengthHistogram();
    readsPointer->writeReadLengthHistogram(fileName);


    cout << "Number of reads: " << readsPointer->readCount() << endl;
    cout << "Number of bases: " << readsPointer->getTotalBaseCount() << endl;
    cout << "Average read length: " <<
        uint64_t(std::round(double(readsPointer->getTotalBaseCount()) / double(readsPointer->readCount()))) << endl;
    cout << "Read N50: " << readsPointer->getN50() << endl;

}



void Assembler::computeReadIdsSortedByName()
{
    readsPointer->computeReadIdsSortedByName();
}
