// Shasta.
#include "Assembler.hpp"
#include "extractKmer.hpp"
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

