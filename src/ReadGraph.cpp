// Shasta.
#include "ReadGraph.hpp"
#include "Anchor.hpp"
#include "deduplicate.hpp"
#include "Reads.hpp"
using namespace shasta2;

// Standard library.
#include "fstream.hpp"

// Explicit instantiation.
#include "MultithreadedObject.tpp"
template class MultithreadedObject<ReadGraph>;


// Initial construction of the ReadGraph from Anchors.
ReadGraph::ReadGraph(
    const Anchors& anchors,
    uint64_t threadCount) :
    MappedMemoryOwner(anchors),
    MultithreadedObject<ReadGraph>(*this),
    anchors(anchors)
{
    // Adjust the numbers of threads, if necessary.
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    const ReadId readCount = anchors.reads.readCount();
    const uint64_t anchorCount = anchors.size();



    // Gather all the OrientedReadIds,
    // including duplicates, encountered for each readId0.

    // Pass 1.
    orientedReadIds.createNew(largeDataName("tmp-ReadGraph-orientedReadIds"), largeDataPageSize);
    orientedReadIds.beginPass1(readCount);
    setupLoadBalancing(anchorCount, 10);
    runThreads(&ReadGraph::threadFunctionPass1, threadCount);

    // Pass 2.
    orientedReadIds.beginPass2();
    setupLoadBalancing(anchorCount, 10);
    runThreads(&ReadGraph::threadFunctionPass2, threadCount);
    orientedReadIds.endPass2(true, true);

    // Pass 3: for each readId0, count the number of times each OrientedReadId appears
    // and generate EdgePairs.
    threadEdgePairs.resize(threadCount);
    setupLoadBalancing(readCount, 10);
    runThreads(&ReadGraph::threadFunctionPass3, threadCount);

    // Gather the EdgePairs found by each thread.
    edgePairs.createNew(largeDataName("ReadGraph-EdgePairs"), largeDataPageSize);
    for(uint64_t threadId=0; threadId<threadCount; threadId++) {
        MemoryMapped::Vector<EdgePair>& thisThreadEdgePairs = *threadEdgePairs[threadId];
        for(const EdgePair& edgePair: thisThreadEdgePairs) {
            edgePairs.push_back(edgePair);
        }
        thisThreadEdgePairs.remove();
        threadEdgePairs[threadId] = 0;
    }

    // Don't keep the orientedReadIds.
    orientedReadIds.remove();

    cout << "The ReadGraph has " << 2 * readCount << " vertices and " <<
        2 * edgePairs.size() << " edges." << endl;

    writeGraphviz();



    // Fill in the connectivity table.
    connectivityTable.createNew(largeDataName("ReadGraph-ConnectivityTable"), largeDataPageSize);
    connectivityTable.beginPass1(readCount);
    setupLoadBalancing(edgePairs.size(), 1000);
    runThreads(&ReadGraph::threadFunctionPass4, threadCount);

    connectivityTable.beginPass2();
    setupLoadBalancing(edgePairs.size(), 1000);
    runThreads(&ReadGraph::threadFunctionPass5, threadCount);
    connectivityTable.endPass2(true, true);
    SHASTA2_ASSERT(connectivityTable.totalSize() == 2 * edgePairs.size());


}



// Access the ReadGraph from binary data.
ReadGraph::ReadGraph(const Anchors& anchors) :
    MappedMemoryOwner(anchors),
    MultithreadedObject<ReadGraph>(*this),
    anchors(anchors)
{
    edgePairs.accessExistingReadOnly(largeDataName("ReadGraph-EdgePairs"));
    connectivityTable.accessExistingReadOnly(largeDataName("ReadGraph-ConnectivityTable"));
}



void ReadGraph::threadFunctionPass1(uint64_t)
{
    threadFunctionPass12(1);
}



void ReadGraph::threadFunctionPass2(uint64_t)
{
    threadFunctionPass12(2);
}



void ReadGraph::threadFunctionPass12(uint64_t pass)
{
    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {
        SHASTA2_ASSERT((begin % 2) == 0);
        SHASTA2_ASSERT((end % 2) == 0);

        // Loop over "positive" (even) AnchorIds assigned to this batch.
        for(AnchorId anchorId=begin; anchorId!=end; anchorId+=2) {
            const Anchor anchor = anchors[anchorId];
            const uint64_t anchorSize = anchor.size();

            // Loop over pairs of AnchorMarkerInfos in this Anchor.
            for(uint64_t i1=1; i1<anchorSize; i1++) {
                const OrientedReadId orientedReadId1 = anchor[i1].orientedReadId;
                const ReadId readId1 = orientedReadId1.getReadId();
                const Strand strand1 = orientedReadId1.getStrand();

                for(uint64_t i0=0; i0<i1; i0++) {
                    const OrientedReadId orientedReadId0 = anchor[i0].orientedReadId;
                    const ReadId readId0 = orientedReadId0.getReadId();
                    const Strand strand0 = orientedReadId0.getStrand();

                    // We are going to store this pair with the lowest numbered ReadId,
                    // readId0, on strand 0.
                    if(pass == 1) {
                        orientedReadIds.incrementCountMultithreaded(readId0);
                    } else {
                        const bool isSameStrand = (strand0 == strand1);
                        const OrientedReadId storedOrientedReadId1(readId1, isSameStrand ? 0 : 1);
                        orientedReadIds.storeMultithreaded(readId0, storedOrientedReadId1);
                    }
                }
            }

        }

    }
}



void ReadGraph::threadFunctionPass3(uint64_t threadId)
{
    // Allocate the EdgePairs found by this thread.
    shared_ptr< MemoryMapped::Vector<EdgePair> >& thisThreadEdgePairsPointer = threadEdgePairs[threadId];
    thisThreadEdgePairsPointer = make_shared< MemoryMapped::Vector<EdgePair> >();
    MemoryMapped::Vector<EdgePair>& thisThreadEdgePairs = *thisThreadEdgePairsPointer;
    thisThreadEdgePairs.createNew(largeDataName("tmp-ReadGraph-EdgePairs-" + to_string(threadId)), largeDataPageSize);

    // Work vectors used below but defined here to reduce
    // memory allocation activity.
    vector<OrientedReadId> orientedReadIds1;
    vector<uint64_t> count;

    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over the ReadIds assigned to this batch.
        for(ReadId readId0=ReadId(begin); readId0!=ReadId(end); readId0++) {

            // Make a copy of the orientedReadIds for this ReadId.
            span<OrientedReadId> orientedReadIds1Span = orientedReadIds[readId0];
            orientedReadIds1.resize(orientedReadIds1Span.size());
            std::ranges::copy(orientedReadIds1Span, orientedReadIds1.begin());

            // Deduplicate and count.
            deduplicateAndCountWithThreshold(orientedReadIds1, count, minCoverage);

            for(uint64_t i=0; i<orientedReadIds1.size(); i++) {
                const OrientedReadId orientedReadId1 = orientedReadIds1[i];
                const bool isSameStrand = orientedReadId1.getStrand() == 0;
                const uint16_t coverage = uint16_t(count[i]);

                EdgePair edgePair(readId0, orientedReadId1.getReadId(), isSameStrand, coverage);
                thisThreadEdgePairs.push_back(edgePair);
            }

        }
    }
}


void ReadGraph::threadFunctionPass4(uint64_t)
{
    threadFunctionPass45(4);
}



void ReadGraph::threadFunctionPass5(uint64_t)
{
    threadFunctionPass45(5);
}



void ReadGraph::threadFunctionPass45(uint64_t pass)
{
    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over all EdgePairs assigned to this batch.
        for(uint64_t edgePairIndex=begin; edgePairIndex!=end; ++edgePairIndex) {
            const EdgePair& edgePair = edgePairs.begin()[edgePairIndex];

            if(pass == 4) {
                connectivityTable.incrementCountMultithreaded(edgePair.readId0);
                connectivityTable.incrementCountMultithreaded(edgePair.readId1);
            } else {
                connectivityTable.storeMultithreaded(edgePair.readId0, edgePairIndex);
                connectivityTable.storeMultithreaded(edgePair.readId1, edgePairIndex);
            }
        }
    }
}



void ReadGraph::writeGraphviz() const
{
    const ReadId readCount = anchors.reads.readCount();

    ofstream dot("ReadGraph.dot");

    dot << "graph ReadGraph {\n";

    // Vertices.
    for(ReadId readId=0; readId<readCount; readId++) {
        for(Strand strand=0; strand<2; strand++) {
            const OrientedReadId orientedReadId(readId, strand);
            dot << "\"" << orientedReadId << "\";\n";
        }
    }

    // Edges.
    for(const EdgePair& edgePair: edgePairs) {

        OrientedReadId orientedReadId0(edgePair.readId0, 0);
        OrientedReadId orientedReadId1(edgePair.readId1, edgePair.isSameStrand ? 0 : 1);
        dot <<  "\"" << orientedReadId0 << "\"--\"" << orientedReadId1 << "\";\n";

        orientedReadId0.flipStrand();
        orientedReadId1.flipStrand();
        dot <<  "\"" << orientedReadId0 << "\"--\"" << orientedReadId1 << "\";\n";

    }

    dot << "}\n";
}
