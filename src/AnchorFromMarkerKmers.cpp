#include "Anchor.hpp"
#include "MarkerKmers.hpp"
#include "performanceLog.hpp"
#include "timestamp.hpp"

using namespace shasta;



Anchors::Anchors(
    const MappedMemoryOwner& mappedMemoryOwner,
    const Reads& reads,
    uint64_t k,
    const Markers& markers,
    shared_ptr<MarkerKmers> markerKmers,
    uint64_t minAnchorCoverage,
    uint64_t maxAnchorCoverage,
    uint64_t threadCount) :
    MultithreadedObject<Anchors>(*this),
    MappedMemoryOwner(mappedMemoryOwner),
    reads(reads),
    k(k),
    markers(markers)
{
    kHalf = k / 2;

    performanceLog << timestamp << "Anchor creation begins." << endl;

    // Store arguments so all threads can see them.
    ConstructFromMarkerKmersData& data = constructFromMarkerKmersData;
    data.minAnchorCoverage = minAnchorCoverage;
    data.maxAnchorCoverage = maxAnchorCoverage;
    data.markerKmers = markerKmers;

    // Use the MarkerKmers to construct anchors.
    // Each thread stores the anchors it finds in a separate vector.
    data.threadAnchors.clear(); // Just in case
    data.threadAnchors.resize(threadCount);
    const uint64_t batchSize = 1000;
    setupLoadBalancing(markerKmers->size(), batchSize);
    runThreads(&Anchors::constructFromMarkerKmersThreadFunction, threadCount);

    // Gather the anchors found by all threads.
    performanceLog << timestamp << "Gathering the anchors found by all threads." << endl;
    anchorMarkerIntervals.createNew(
            largeDataName("AnchorMarkerIntervals"),
            largeDataPageSize);
    anchorSequences.createNew(
        largeDataName("AnchorSequences"), largeDataPageSize);
    anchorInfos.createNew(largeDataName("AnchorInfos"), largeDataPageSize);
    for(uint64_t threadId=0; threadId<threadCount; threadId++) {
        auto& threadAnchorsPointer = data.threadAnchors[threadId];
        auto& threadAnchors = *threadAnchorsPointer;

        // Loop over the anchors found by this thread.
        for(uint64_t i=0; i<threadAnchors.size(); i++) {
            const auto threadAnchor = threadAnchors[i];
            anchorMarkerIntervals.appendVector();
            for(const auto& markerInfo: threadAnchor) {
                anchorMarkerIntervals.append(AnchorMarkerInterval(markerInfo.orientedReadId, markerInfo.ordinal));
            }
        }
        threadAnchors.remove();
        threadAnchorsPointer = 0;
    }

    // The anchor sequences are all empty.
    // The ordinal offsets are all 0.
    const uint64_t anchorCount = anchorMarkerIntervals.size();
    anchorInfos.resize(anchorCount);
    for(AnchorId anchorId=0; anchorId<anchorCount; anchorId++) {
        anchorSequences.appendVector();
        anchorInfos[anchorId].ordinalOffset = 0;
    }

    cout << "Number of anchors per strand: " << anchorCount / 2 << endl;
    performanceLog << timestamp << "Anchor creation from marker kmers ends." << endl;
}



void Anchors::constructFromMarkerKmersThreadFunction(uint64_t threadId)
{
    using MarkerInfo = MarkerKmers::MarkerInfo;

    ConstructFromMarkerKmersData& data = constructFromMarkerKmersData;
    const uint64_t minAnchorCoverage = data.minAnchorCoverage;
    const uint64_t maxAnchorCoverage = data.maxAnchorCoverage;
    const MarkerKmers& markerKmers = *data.markerKmers;

    // Initialize the anchors that will be found by this thread.
    auto& threadAnchorsPointer = data.threadAnchors[threadId];
    threadAnchorsPointer = make_shared<MemoryMapped::VectorOfVectors<ConstructFromMarkerKmersData::MarkerInfo, uint64_t> >();
    auto& threadAnchors = *threadAnchorsPointer;
    threadAnchors.createNew(largeDataName("tmp-threadAnchors-" + to_string(threadId)), largeDataPageSize);

    // A vector used below and defined here to reduce memory allocation activity.
    // It will contain the MarkerInfos for a marker Kmer, ecluding
    // the ones for which the same ReadId appears more than once in the same Kmer.
    // There are the ones that will be used to generate anchors.
    vector<MarkerInfo> usableMarkerInfos;

    // Loop over batches of marker Kmers assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over marker Kmers assigned to this batch.
        for(uint64_t i=begin; i!=end; i++) {

            // Get the MarkerInfos for this marker Kmer.
            const span<const MarkerInfo> markerInfos = markerKmers[i];

            // Check for high coverage using all of the marker infos.
            if(markerInfos.size() > maxAnchorCoverage) {
                continue;
            }



            // Gather the usable MarkerInfos.
            // These are the ones for which the ReadId is different from the ReadId
            // of the previous and next MarkerInfo.
            usableMarkerInfos.clear();
            for(uint64_t i=0; i<markerInfos.size(); i++) {
                const MarkerInfo& markerInfo = markerInfos[i];
                bool isUsable = true;

                // Check if same ReadId of previous MarkerInfo.
                if(i != 0) {
                    isUsable =
                        isUsable and
                        (markerInfo.orientedReadId.getReadId() != markerInfos[i-1].orientedReadId.getReadId());
                }

                // Check if same ReadId of next MarkerInfo.
                if(i != markerInfos.size() - 1) {
                    isUsable =
                        isUsable and
                        (markerInfo.orientedReadId.getReadId() != markerInfos[i+1].orientedReadId.getReadId());
                }

                if(isUsable) {
                    usableMarkerInfos.push_back(markerInfo);
                }
            }



            // Check for low coverage using the usable marker infos.
            if(usableMarkerInfos.size() < minAnchorCoverage) {
                continue;
            }



            // If getting here, we will generate a pair of Anchors corresponding
            // to this Kmer.


            // Generate the first anchor of the pair (no reverse complementing).
            threadAnchors.appendVector(usableMarkerInfos);


            // Reverse complement the usableMarkerInfos, then
            // generate the second anchor in the pair.
            for(MarkerInfo& markerInfo: usableMarkerInfos) {
                markerInfo = markerKmers.reverseComplement(markerInfo);
            }
            threadAnchors.appendVector(usableMarkerInfos);
        }
    }

}

