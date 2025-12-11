#include "ReadGraph.hpp"
#include "Anchor.hpp"
#include "deduplicate.hpp"
#include "Reads.hpp"
using namespace shasta2;

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

    // Pass3: for each readId0, count the number of times each OrientedReadId appears,
    // and figure out how many edges we will generate for that readId0.
    edges.createNew(largeDataName("ReadGraph-Edges"), largeDataPageSize);
    edges.beginPass1(readCount);
    setupLoadBalancing(readCount, 10);
    runThreads(&ReadGraph::threadFunctionPass3, threadCount);

    // Pass4: do it again, this time storing the edges.
    edges.beginPass2();
    setupLoadBalancing(readCount, 10);
    runThreads(&ReadGraph::threadFunctionPass4, threadCount);
    edges.endPass2(true, true);

    // Don't keep the orientedReadIds.
    orientedReadIds.remove();

    cout << "The ReadGraph has " << 2 * readCount << " vertices and " <<
        2 * edges.totalSize() << " edges." << endl;
}



// Access the ReadGraph from binary data.
ReadGraph::ReadGraph(const Anchors& anchors) :
    MappedMemoryOwner(anchors),
    MultithreadedObject<ReadGraph>(*this),
    anchors(anchors)
{

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



void ReadGraph::threadFunctionPass3(uint64_t)
{
    threadFunctionPass34(3);
}



void ReadGraph::threadFunctionPass4(uint64_t)
{
    threadFunctionPass34(4);
}



void ReadGraph::threadFunctionPass34(uint64_t pass)
{
    // Work vector used below but defined here to reduce
    // memory allocation activity.
    vector<OrientedReadId> v;
    vector<uint64_t> count;

    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over the ReadIds assigned to this batch.
        for(ReadId readId=ReadId(begin); readId!=ReadId(end); readId++) {

            // Make a copy of the orientedReadIds for this ReadId.
            span<OrientedReadId> s = orientedReadIds[readId];
            v.resize(s.size());
            std::ranges::copy(s, v.begin());

            // Deduplicate and count.
            deduplicateAndCountWithThreshold(v, count, minCoverage);

            if(pass == 3) {
                edges.incrementCountMultithreaded(readId, ReadId(count.size()));
            } else {
                for(uint64_t i=0; i<count.size(); i++) {
                    Edge edge(v[i], uint32_t(count[i]));
                    edges.storeMultithreaded(readId, edge);
                }
            }

        }
    }
}
