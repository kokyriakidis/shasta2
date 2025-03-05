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


    cout << "Read statistics for reads that will be used in this assembly:" << endl;
    cout << "    Total number of reads is " << readsPointer->readCount() << "." << endl;
    cout << "    Total number of bases is " << readsPointer->getTotalBaseCount() << "." << endl;
    cout << "    Average read length is " << double(readsPointer->getTotalBaseCount()) / double(readsPointer->readCount());
    cout << " bases." << endl;
    cout << "    N50 for read length is " << readsPointer->getN50() << " bases." << endl;

}



void Assembler::computeReadIdsSortedByName()
{
    readsPointer->computeReadIdsSortedByName();
}
