#include "Journeys.hpp"
#include "Markers.hpp"
#include "orderPairs.hpp"
#include "performanceLog.hpp"
#include "ReadId.hpp"
#include "Reads.hpp"
#include "timestamp.hpp"
using namespace shasta;

// Explicit instantiation.
#include "MultithreadedObject.tpp"
template class MultithreadedObject<Journeys>;



// Initial creation.
// This sets the positionInJourney for every AnchorMarkerInfo
// stored in the Anchors, and for this reason the Anchors
// are not passed in as const.
Journeys::Journeys(
    uint64_t orientedReadCount,
    shared_ptr<Anchors> anchorsPointer,
    uint64_t threadCount,
    const MappedMemoryOwner& mappedMemoryOwner) :
    MultithreadedObject<Journeys>(*this),
    MappedMemoryOwner(mappedMemoryOwner),
    anchorsPointer(anchorsPointer)

{
    performanceLog << timestamp << "Journeys creation begins." << endl;

    // Adjust the numbers of threads, if necessary.
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    const uint64_t anchorCount = anchorsPointer->size();
    const uint64_t anchorBatchCount = 1000;
    const uint64_t orientedReadBatchCount = 1000;

    // Pass1: make space for the journeysWithOrdinals.
    journeysWithOrdinals.createNew(largeDataName("tmp-JourneysWithOrdinals"), largeDataPageSize);
    journeysWithOrdinals.beginPass1(orientedReadCount);
    setupLoadBalancing(anchorCount, anchorBatchCount);
    runThreads(&Journeys::threadFunction1, threadCount);

    // Pass2: store the unsorted journeysWithOrdinals.
    journeysWithOrdinals.beginPass2();
    setupLoadBalancing(anchorCount, anchorBatchCount);
    runThreads(&Journeys::threadFunction2, threadCount);
    journeysWithOrdinals.endPass2();

    // Pass 3:sort the journeysWithOrdinals and make space for the journeys
    journeys.createNew(largeDataName("Journeys"), largeDataPageSize);
    journeys.beginPass1(orientedReadCount);
    setupLoadBalancing(orientedReadCount, orientedReadBatchCount);
    runThreads(&Journeys::threadFunction3, threadCount);

    // Pass 4: copy the sorted journeysWithOrdinals to the journeys.
    journeys.beginPass2();
    setupLoadBalancing(orientedReadCount, orientedReadBatchCount);
    runThreads(&Journeys::threadFunction4, threadCount);
    journeys.endPass2(false, true);

    journeysWithOrdinals.remove();

    performanceLog << timestamp << "Journeys creation ends." << endl;
}



void Journeys::threadFunction1(uint64_t /* threadId */)
{
    threadFunction12(1);
}



void Journeys::threadFunction2(uint64_t /* threadId */)
{
    threadFunction12(2);
}



void Journeys::threadFunction12(uint64_t pass)
{
    const Anchors& anchors = *anchorsPointer;

    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over all AnchorIds in this batch.
        for(AnchorId anchorId=begin; anchorId!=end; anchorId++) {
            Anchor anchor = anchors[anchorId];

            // Loop over the marker intervals of this Anchor.
            for(const auto& anchorMarkerInterval: anchor) {
                const auto orientedReadIdValue = anchorMarkerInterval.orientedReadId.getValue();

                if(pass == 1) {
                    journeysWithOrdinals.incrementCountMultithreaded(orientedReadIdValue);
                } else {
                    journeysWithOrdinals.storeMultithreaded(
                        orientedReadIdValue, {anchorId, anchorMarkerInterval.ordinal});
                }
            }
        }
    }
}



void Journeys::threadFunction3(uint64_t /* threadId */)
{
    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over all oriented reads assigned to this thread.
        for(uint64_t orientedReadValue=begin; orientedReadValue!=end; orientedReadValue++) {
            auto v = journeysWithOrdinals[orientedReadValue];
            sort(v.begin(), v.end(), OrderPairsBySecondOnly<uint64_t, uint32_t>());
            journeys.incrementCountMultithreaded(orientedReadValue, v.size());
        }
    }
}



void Journeys::threadFunction4(uint64_t /* threadId */)
{
    Anchors& anchors = *anchorsPointer;

    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over all oriented reads assigned to this thread.
        for(uint64_t orientedReadValue=begin; orientedReadValue!=end; orientedReadValue++) {
            const OrientedReadId orientedReadId = OrientedReadId::fromValue(ReadId(orientedReadValue));

            // Copy the journeysWithOrdinals to the journeys.
            const auto v = journeysWithOrdinals[orientedReadValue];
            const auto journey = journeys[orientedReadValue];
            SHASTA_ASSERT(journey.size() == v.size());
            for(uint64_t i=0; i<v.size(); i++) {
                journey[i] = v[i].first;
            }

            // Store journey information for this oriented read in the marker interval.
            for(uint64_t position=0; position<journey.size(); position++) {
                const AnchorId anchorId = journey[position];
                span<AnchorMarkerInfo> markerInfos = anchors.anchorMarkerInfos[anchorId];
                bool found = false;
                for(AnchorMarkerInfo& markerInfo: markerInfos) {
                    if(markerInfo.orientedReadId == orientedReadId) {
                        markerInfo.positionInJourney = uint32_t(position);
                        found = true;
                        break;
                    }
                }
                SHASTA_ASSERT(found);
            }
        }
    }
}



// For each read, write out the largest gap between adjacent anchors.
// The two oriented reads for a read have the same gaps.
void Journeys::writeAnchorGapsByRead(
    const Reads& reads,
    const Markers& markers,
    const Anchors& anchors
    ) const
{

    // Open the output csv file and write a header.
    ofstream csv("AnchorGaps.csv");
    csv << "ReadId,Length,Gap\n";

    // Loop over all reads.
    for(ReadId readId=0; readId<reads.readCount(); readId++) {
        const uint64_t readLength = reads.getRead(readId).baseCount;

        // Put in on strand 0. It would have the same gap on strand 1.
        const OrientedReadId orientedReadId(readId, 0);

        // Get the markers and the journey.
        const auto orientedReadMarkers = markers[orientedReadId.getValue()];
        const auto journey = (*this)[orientedReadId];

        // Loop over adjacent positions in journey.
        uint64_t previousPosition = 0;
        uint64_t maxGap = 0;
        for(uint64_t i=0; i<=journey.size(); i++) {

            uint64_t position = invalid<uint64_t>;
            if(i == journey.size()) {
                position = readLength;
            } else {
                const AnchorId anchorId = journey[i];
                const Anchor anchor = anchors[anchorId];
                for(const AnchorMarkerInfo& markerInfo: anchor) {
                    if(markerInfo.orientedReadId == orientedReadId) {
                        const uint32_t ordinal = markerInfo.ordinal;
                        const Marker& marker = orientedReadMarkers[ordinal];
                        position = marker.position;
                        break;
                    }
                }
            }
            SHASTA_ASSERT(position != invalid<uint64_t>);
            SHASTA_ASSERT(position >= previousPosition);

            const uint64_t gap = position - previousPosition;
            maxGap = max(maxGap, gap);

            previousPosition = position;
        }
        csv << readId << "," << readLength << "," << maxGap << "\n";
    }
}



// Access from binary data.
Journeys::Journeys(const MappedMemoryOwner& mappedMemoryOwner) :
    MultithreadedObject<Journeys>(*this),
    MappedMemoryOwner(mappedMemoryOwner)
{
    journeys.accessExistingReadOnly(largeDataName("Journeys"));
}
