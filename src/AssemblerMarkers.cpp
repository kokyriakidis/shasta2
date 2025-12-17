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
    const ReadId readCount = reads().readCount();

    ofstream csv("PalindromicMetrics.csv");
    csv << "ReadId,PalindromicRate\n";

    // Loop over all reads.
    for(ReadId readId=0; readId<readCount; readId++) {

        const double palindromicRate = analyzeStrandReversal(readId, false);

        csv << readId << ",";
        csv << palindromicRate << ",";
        csv << "\n";
    }
}



double Assembler::analyzeStrandReversal(ReadId readId, bool debug) const
{
    // EXPOSE WHEN CODE STABILIZES.
    const double driftRate = 0.01;

    const uint64_t k = assemblerInfo->k;

    // Access the sequence of this read (without reverse complementing).
    const LongBaseSequenceView sequence0 = reads().getRead(readId);
    const uint32_t length = uint32_t(sequence0.baseCount);
    const uint32_t maxDrift = uint32_t(std::round(driftRate * double(length)));

    if(debug) {
        cout << "Length " << length << endl;
        cout << "Max drift " << maxDrift << endl;
    }

    // Access the markers of this read (without reverse complementing.
    const auto markers0 = markers()[OrientedReadId(readId, 0).getValue()];


    // Gather all marker Kmers of this read and its reverse complement,
    // and the positions in which each Kmer appears.
    std::map<Kmer, vector<uint32_t> > kmerTable0;
    std::map<Kmer, vector<uint32_t> > kmerTable1;

    for(const Marker& marker0: markers0) {
        const uint32_t position0 = marker0.position;
        Kmer kmer0;
        extractKmer128(sequence0, position0, k, kmer0);
        kmerTable0[kmer0].push_back(position0);

        const uint32_t position1 = length - 1 - position0;
        const Kmer kmer1 = kmer0.reverseComplement(k);
        kmerTable1[kmer1].push_back(position1);
    }

    for(auto& p: kmerTable1) {
        vector<uint32_t>& positions1 = p.second;
        std::ranges::reverse(positions1);
    }


    if(debug) {
        // Loop over common marker k-mers.
        ofstream csv("AnalyzeStrandReversal.csv");
        for(const auto& p0: kmerTable0) {
            const Kmer& kmer = p0.first;
            auto it1 = kmerTable1.find(kmer);
            if(it1 != kmerTable1.end()) {
                const auto& p1 = *it1;
                const vector<uint32_t>& positions0 = p0.second;
                const vector<uint32_t>& positions1 = p1.second;
                for(const uint32_t position0: positions0) {
                    for(const uint32_t position1: positions1) {
                        csv <<
                            position0 << "," <<
                            position1 << "," <<
                            int32_t(position1) - int32_t(position0) << "\n";
                    }
                }
            }
        }
    }


    // To decide if this read is palindromic, loop over all its marker kmers.
    // Find how many times we find the same k-mer at a nearby position
    // in the reverse complemented read.
    uint64_t totalCount = 0;
    uint64_t successCount = 0;
    for(const auto& p0: kmerTable0) {
        const Kmer& kmer = p0.first;
        const vector<uint32_t>& positions0 = p0.second;
        totalCount += positions0.size();

        const auto it1 = kmerTable1.find(kmer);
        if(it1 != kmerTable1.end()) {
            const vector<uint32_t>& positions1 = it1->second;
            SHASTA2_ASSERT(not positions1.empty());

            // This can be made more efficient.
            for(uint64_t position0: positions0) {
                for(uint64_t position1: positions1) {
                    if(abs(int32_t(position0) - int32_t(position1)) < maxDrift) {
                        ++successCount;
                        break;
                    }
                }
            }
        }
    }

    const double successRate = double(successCount) / double(totalCount);

    if(debug) {
        cout << "Total " << totalCount << ", success " << successCount <<
            ", success rate " << successRate << endl;
    }

    return successRate;

}
