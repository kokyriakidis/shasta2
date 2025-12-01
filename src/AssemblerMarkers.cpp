// Shasta.
#include "Assembler.hpp"
#include "deduplicate.hpp"
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



void Assembler::findPalindromicReads() const
{
    const uint64_t k = assemblerInfo->k;
    const ReadId readCount = reads().readCount();

    ofstream csv("PalindromicMetrics.csv");
    csv << "ReadId,Intersection,Union,Jaccard,\n";

    // Work vectors used below but defined here to reduce
    // memory allocation activity.
    vector<Kmer> kmers0;
    vector<Kmer> kmers1;
    vector<uint64_t> count0;
    vector<uint64_t> count1;

    // Loop over all reads.
    for(ReadId readId=0; readId<readCount; readId++) {

        // Access the sequence of this read (without reverse complementing).
        const LongBaseSequenceView sequence0 = reads().getRead(readId);

        // Access the markers of this read (without reverse complementing.
        const auto markers0 = markers()[OrientedReadId(readId, 0).getValue()];

        // Gather all marker Kmers of this read and its reverse complement.
        kmers0.clear();
        kmers1.clear();
        for(const Marker& marker0: markers0) {
            const uint32_t position = marker0.position;
            Kmer kmer0;
            extractKmer128(sequence0, position, k, kmer0);
            kmers0.push_back(kmer0);
            const Kmer kmer1 = kmer0.reverseComplement(k);
            kmers1.push_back(kmer1);
        }

        // Deduplicate and count. Effectively we compute multisets (bags) of Kmers.
        deduplicateAndCount(kmers0, count0);
        deduplicateAndCount(kmers1, count1);



        // Compute the Jaccard similarity of the two multisets (bags).
        // See https://en.wikipedia.org/wiki/Multiset
        // and https://en.wikipedia.org/wiki/Jaccard_index
        // (in the second link, look for "multisets" under "Overview".

        // For the union we just need the total size of the two multisets.
        const uint64_t sum0 = std::accumulate(count0.begin(), count0.end(), 0);
        const uint64_t sum1 = std::accumulate(count1.begin(), count1.end(), 0);
        const uint64_t unionSize = sum0 + sum1;

        // For the intersection we do a joint loop over the common Kmers,
        // which are now sorted.
        uint64_t intersectionSize = 0;
        uint64_t i0 = 0;
        uint64_t i1 = 0;
        while((i0 < kmers0.size()) and (i1 < kmers1.size())) {
            const Kmer& kmer0 = kmers0[i0];
            const Kmer& kmer1 = kmers1[i1];
            if(kmer0 < kmer1) {
                ++i0;
            } else if(kmer1 < kmer0) {
                ++i1;
            } else {
                intersectionSize += min(count0[i0], count1[i1]);
                ++i0;
                ++i1;
            }
        }

        const double jaccard = double(intersectionSize) / double(unionSize);
        csv << readId << ",";
        csv << intersectionSize << ",";
        csv << unionSize << ",";
        csv << jaccard << ",";
        csv << "\n";
    }
}
