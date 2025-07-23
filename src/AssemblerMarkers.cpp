// Shasta.
#include "Assembler.hpp"
#include "extractKmer128.hpp"
#include "Markers.hpp"
#include "MarkerKmers.hpp"
#include "performanceLog.hpp"
#include "Reads.hpp"
#include "timestamp.hpp"
using namespace shasta;

// Standard library.
#include "fstream.hpp"


void Assembler::createMarkers(size_t threadCount)
{
    readsPointer->checkReadsAreOpen();
    SHASTA_ASSERT(kmerChecker);

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



void Assembler::createMarkerKmers(uint64_t threadCount)
{
    const MappedMemoryOwner& mappedMemoryOwner = *this;

    markerKmers = make_shared<MarkerKmers>(
        assemblerInfo->k,
        mappedMemoryOwner,
        reads(),
        markers(),
        threadCount);

    computeMarkerErrorRates();
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
void Assembler::computeMarkerErrorRates() const
{
    SHASTA_ASSERT(markerKmers);
    const uint64_t k = assemblerInfo->k;

    ofstream csv("MarkerErrorRates.csv");
    csv << "ReadId,Length,TotalMarkerCount,FrequencyOneMarkerCount,FrequencyOneMarkerRatio,\n";

    // Loop over all reads.
    for(ReadId readId=0; readId<reads().readCount(); readId++) {
        const auto readSequence = reads().getRead(readId);
        const auto readMarkers = markers()[OrientedReadId(readId, 0).getValue()];

        // Count frequency 1 marker kmers in this read.
        uint64_t count = 0;
        for(const Marker& marker: readMarkers) {
            const uint32_t position = marker.position;
            Kmer kmer;
            extractKmer128(readSequence, position, k, kmer);
            const uint64_t frequency = markerKmers->getFrequency(kmer);
            if(frequency == 1) {
                ++count;
            }
        }
        csv << readId << ",";
        csv << readSequence.baseCount << ",";
        csv << readMarkers.size() << ",";
        csv << count << ",";
        csv << double(count) / double(readMarkers.size()) << ",";
        csv << "\n";

    }
}


