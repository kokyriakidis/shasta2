// Shasta.
#include "Assembler.hpp"
#include "extractKmer128.hpp"
#include "Markers.hpp"
#include "MarkerKmers.hpp"
#include "performanceLog.hpp"
#include "Reads.hpp"
#include "ReadSummary.hpp"
#include "timestamp.hpp"
using namespace shasta2;

// Standard library.
#include "fstream.hpp"


void Assembler::createMarkers(size_t threadCount)
{
    readsPointer->checkReadsAreOpen();
    SHASTA2_ASSERT(kmerChecker);

    markersPointer = make_shared<Markers>(
        *this,
        assemblerInfo->k,
        readsPointer,
        kmerChecker,
        threadCount);
}



void Assembler::accessMarkers()
{
    markersPointer = make_shared<Markers>(*this, assemblerInfo->k, readsPointer);
}



void Assembler::checkMarkersAreOpen() const
{
    if(not(markersPointer and !markersPointer->isOpen())) {
        throw runtime_error("Markers are not accessible.");
    }
}



void Assembler::createMarkerKmers(double maxMarkerErrorRate, uint64_t threadCount)
{

    const MappedMemoryOwner& mappedMemoryOwner = *this;
    const uint64_t readCount = reads().readCount();

    // First, create marker k-mers using all reads.
    vector<bool> useRead(readCount, true);
    markerKmers = make_shared<MarkerKmers>(
        assemblerInfo->k,
        mappedMemoryOwner,
        reads(),
        useRead,
        markers(),
        threadCount);

    // Compute marker error rates using all reads and store them in the ReadSummaries.
    // Also set the isUsedForAssembly flag.
    vector<uint64_t> lowFrequencyMarkerCount;
    computeMarkerErrorRates(lowFrequencyMarkerCount);
    for(ReadId readId=0; readId<readCount; readId++) {
        const auto readMarkers = markers()[OrientedReadId(readId, 0).getValue()];
        const uint64_t readMarkerCount = readMarkers.size();
        const double markerErrorRate = double(lowFrequencyMarkerCount[readId]) / double(readMarkerCount);
        useRead[readId] = (markerErrorRate <= maxMarkerErrorRate);
        ReadSummary& readSummary = readSummaries[readId];
        readSummary.isUsedForAssembly = useRead[readId];
        readSummary.initialMarkerErrorRate = markerErrorRate;

    }


    // Count the reads and bases that will be used for assembly.
    uint64_t keepCount = 0;
    uint64_t baseKeepCount = 0;
    for(ReadId readId=0; readId<readCount; readId++) {
        if(useRead[readId]) {
            ++keepCount;
            baseKeepCount += reads().getRead(readId).baseCount;
        }
    }
    cout << "Kept " << keepCount << " reads with low marker error rate out of " <<
        readCount << " total (" << baseKeepCount <<
        " bases out of " << reads().getTotalBaseCount() << " total)." << endl;

    // Do it again, this time using only the reads with low marker error rate.
    markerKmers->remove();
    markerKmers = make_shared<MarkerKmers>(
        assemblerInfo->k,
        mappedMemoryOwner,
        reads(),
        useRead,
        markers(),
        threadCount);

    // Compute marker error rates.
    computeMarkerErrorRates(lowFrequencyMarkerCount);
    for(ReadId readId=0; readId<readCount; readId++) {
        const auto readMarkers = markers()[OrientedReadId(readId, 0).getValue()];
        const uint64_t readMarkerCount = readMarkers.size();
        const double markerErrorRate = double(lowFrequencyMarkerCount[readId]) / double(readMarkerCount);
        ReadSummary& readSummary = readSummaries[readId];
        readSummary.markerErrorRate = markerErrorRate;
    }
 }



void Assembler::accessMarkerKmers()
{
    const MappedMemoryOwner& mappedMemoryOwner = *this;

    markerKmers = make_shared<MarkerKmers>(
        assemblerInfo->k,
        mappedMemoryOwner,
        reads(),
        markers());
}



// Compute marker error rates for each read by counting the
// marker k-mers with frequency 1.
void Assembler::computeMarkerErrorRates(
    vector<uint64_t>& lowFrequencyMarkerCount) const
{
    SHASTA2_ASSERT(markerKmers);
    const uint64_t k = assemblerInfo->k;
    const uint64_t readCount = reads().readCount();

    lowFrequencyMarkerCount.clear();
    lowFrequencyMarkerCount.resize(readCount, 0);

    // Loop over all reads.
    for(ReadId readId=0; readId<readCount; readId++) {
        const auto readSequence = reads().getRead(readId);
        const auto readMarkers = markers()[OrientedReadId(readId, 0).getValue()];

        // Count frequency 1 marker kmers in this read.
        for(const Marker& marker: readMarkers) {
            const uint32_t position = marker.position;
            Kmer kmer;
            extractKmer128(readSequence, position, k, kmer);
            const uint64_t frequency = markerKmers->getFrequency(kmer);
            // Also count frequency 0.
            // Possible, for reads for which useRead is false;
            if(frequency <= 1) {
                ++lowFrequencyMarkerCount[readId];
            }
        }
    }
}
