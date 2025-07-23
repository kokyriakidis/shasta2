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



void Assembler::createMarkerKmers(double maxMarkerErrorRate, uint64_t threadCount)
{

    const MappedMemoryOwner& mappedMemoryOwner = *this;
    const uint64_t readCount = reads().readCount();

    // First, compute marker k-mers using all reads.
    vector<bool> useRead(readCount, true);
    markerKmers = make_shared<MarkerKmers>(
        assemblerInfo->k,
        mappedMemoryOwner,
        reads(),
        useRead,
        markers(),
        threadCount);

    // Compute marker error rates.
    vector<uint64_t> frequencyOneMarkerCount;
    computeMarkerErrorRates(useRead, frequencyOneMarkerCount);
    writeMarkerErrorRates(useRead, frequencyOneMarkerCount, "MarkerErrorRates-AllReads.csv");



    // Do it again, this time using only the reads with low marker error rate.
    for(ReadId readId=0; readId<readCount; readId++) {
        const auto readMarkers = markers()[OrientedReadId(readId, 0).getValue()];
        const double markerErrorRate = double(frequencyOneMarkerCount[readId]) / double(readMarkers.size());
        useRead[readId] = (markerErrorRate <= maxMarkerErrorRate);
    }

    const uint64_t keepCount = std::ranges::count(useRead, 1);
    uint64_t baseKeepCount = 0;
    for(ReadId readId=0; readId<readCount; readId++) {
        if(useRead[readId]) {
            baseKeepCount += reads().getRead(readId).baseCount;
        }
    }
    cout << "Kept " << keepCount << " reads with low marker error rate out of " <<
        readCount << " total (" << baseKeepCount <<
        " bases out of " << reads().getTotalBaseCount() << " total)." << endl;

    markerKmers->remove();
    markerKmers = make_shared<MarkerKmers>(
        assemblerInfo->k,
        mappedMemoryOwner,
        reads(),
        useRead,
        markers(),
        threadCount);

    // Compute marker error rates.
    computeMarkerErrorRates(useRead, frequencyOneMarkerCount);
    writeMarkerErrorRates(useRead, frequencyOneMarkerCount, "MarkerErrorRates.csv");
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
    const vector<bool>& useRead,
    vector<uint64_t>& frequencyOneMarkerCount) const
{
    SHASTA_ASSERT(markerKmers);
    const uint64_t k = assemblerInfo->k;
    const uint64_t readCount = reads().readCount();

    frequencyOneMarkerCount.clear();
    frequencyOneMarkerCount.resize(readCount, 0);

    // Loop over all reads.
    for(ReadId readId=0; readId<readCount; readId++) {
        if(not useRead[readId]) {
            continue;
        }
        const auto readSequence = reads().getRead(readId);
        const auto readMarkers = markers()[OrientedReadId(readId, 0).getValue()];

        // Count frequency 1 marker kmers in this read.
        for(const Marker& marker: readMarkers) {
            const uint32_t position = marker.position;
            Kmer kmer;
            extractKmer128(readSequence, position, k, kmer);
            const uint64_t frequency = markerKmers->getFrequency(kmer);
            if(frequency == 1) {
                ++frequencyOneMarkerCount[readId];
            }
        }
    }
}



void Assembler::writeMarkerErrorRates(
    const vector<bool>& useRead,
    const vector<uint64_t>& frequencyOneMarkerCount,
    const string& fileName) const
{
    ofstream csv(fileName);
    csv << "ReadId,Length,TotalMarkerCount,FrequencyOneMarkerCount,FrequencyOneMarkerRatio,\n";

    // Loop over all reads.
    for(ReadId readId=0; readId<reads().readCount(); readId++) {
        if(not useRead[readId]) {
            continue;
        }
        const auto readSequence = reads().getRead(readId);
        const auto readMarkers = markers()[OrientedReadId(readId, 0).getValue()];

        csv << readId << ",";
        csv << readSequence.baseCount << ",";
        csv << readMarkers.size() << ",";
        csv << frequencyOneMarkerCount[readId] << ",";
        csv << double(frequencyOneMarkerCount[readId]) / double(readMarkers.size()) << ",";
        csv << "\n";

    }
}
