// Shasta.
#include "ReadGraph.hpp"
#include "Anchor.hpp"
#include "deduplicate.hpp"
#include "DisjointSets.hpp"
#include "Reads.hpp"
using namespace shasta2;

// Standard library.
#include "fstream.hpp"
#include <queue>
#include <random>
#include "tuple.hpp"

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
    cout << "The ReadGraph has " << 2 * readCount << " vertices and " <<
        2 * edgePairs.size() << " edges." << endl;

    // Don't keep the orientedReadIds.
    orientedReadIds.remove();


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

    setupLoadBalancing(readCount, 100);
    runThreads(&ReadGraph::threadFunctionPass6, threadCount);

    findStrandSymmetriQuadrilaterals();

#if 0
    // Use self-complementary paths flag cross-strand EdgePairs.
    for(uint64_t maxDistance=2; maxDistance<=2; maxDistance++) {
        for(uint64_t iteration=0; iteration<1; iteration++) {
            if(flagCrossStrandEdgePairs(maxDistance) ==0) {
                break;
            }
        }
    }
#endif

#if 0
    // Compute the connected components of the complete, double-stranded ReadGraph,
    // excluding EdgePairs flagged as cross-strand
    computeComponents2();

    // Compute the connected components of a partial, single-stranded ReadGraph.
    computeComponents1();
#endif

#if 0
    writeConnectivityTable();
    writeGraphviz();
#endif

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

                // Put the lowest ReadId as the first in the edgePair.
                ReadId readIdA = readId0;
                ReadId readIdB = orientedReadId1.getReadId();
                if(readIdB < readIdA) {
                    swap(readIdA, readIdA);
                }

                EdgePair edgePair(readIdA, readIdB, isSameStrand, coverage);
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



// This sort the connectivityTable entries for each read.
void ReadGraph::threadFunctionPass6(uint64_t)
{
    // Class used below to sort connectivity table entries.
    class Info {
    public:
        uint64_t edgePairIndex;
        ReadId otherReadId;
        bool isSameStrand;

        Info(
            uint64_t edgePairIndex,
            ReadId otherReadId,
            bool isSameStrand) :
            edgePairIndex(edgePairIndex),
            otherReadId(otherReadId),
            isSameStrand(isSameStrand)
        {}

        bool operator<(const Info& that) const
        {
            return tie(otherReadId, isSameStrand) < tie(that.otherReadId, that.isSameStrand);
        }
    };
    vector<Info> infos;



    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over all ReadIds assigned to this batch.
        for(uint64_t readId=begin; readId!=end; ++readId) {

            // Gather the Infos for this ReadId.
            span<uint64_t> edgePairIndexes = connectivityTable[readId];
            infos.clear();
            for(uint64_t edgePairIndex: edgePairIndexes) {
                const EdgePair& edgePair = edgePairs[edgePairIndex];
                infos.emplace_back(edgePairIndex, edgePair.getOther(ReadId(readId)), edgePair.isSameStrand);
            }

            // Sort them.
            sort(infos.begin(), infos.end());

            // Overwrite the edgePairIndexes.
            for(uint64_t i=0; i<infos.size(); i++) {
                const Info& info = infos[i];
                edgePairIndexes[i] = info.edgePairIndex;
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
            const uint64_t componentId =
                (components1.componentId.empty() ? invalid<uint64_t> :
                components1.componentId[orientedReadId.getValue()]);

            dot << "\"" << orientedReadId << "\"";
            if(componentId != invalid<uint64_t>) {
                dot << " [color=green]";
            }
            dot << ";\n";
        }
    }

    // Edges.
    for(const EdgePair& edgePair: edgePairs) {
        if(edgePair.isCrossStrand) {
            continue;
        }

        OrientedReadId orientedReadId0(edgePair.readId0, 0);
        OrientedReadId orientedReadId1(edgePair.readId1, edgePair.isSameStrand ? 0 : 1);
        dot <<  "\"" << orientedReadId0 << "\"--\"" << orientedReadId1 << "\"";
        if(edgePair.isCrossStrand) {
            dot << " [color=yellow]";
        }
        dot << ";\n";

        orientedReadId0.flipStrand();
        orientedReadId1.flipStrand();
        dot <<  "\"" << orientedReadId0 << "\"--\"" << orientedReadId1 << "\"";
        if(edgePair.isCrossStrand) {
            dot << " [color=yellow]";
        }
        dot << ";\n";

    }

    dot << "}\n";
}



void ReadGraph::computeComponents2()
{
    const ReadId readCount = anchors.reads.readCount();
    const uint64_t orientedReadCount = 2 * readCount;

    DisjointSets disjointSets(orientedReadCount);
    for(const EdgePair& edgePair: edgePairs) {
        if(edgePair.isCrossStrand) {
            continue;
        }
        OrientedReadId orientedReadId0(edgePair.readId0, 0);
        OrientedReadId orientedReadId1(edgePair.readId1, edgePair.isSameStrand ? 0 : 1);

        disjointSets.unionSet(orientedReadId0.getValue(), orientedReadId1.getValue());

        orientedReadId0.flipStrand();
        orientedReadId1.flipStrand();

        disjointSets.unionSet(orientedReadId0.getValue(), orientedReadId1.getValue());
    }

    vector< vector<uint64_t> > rawComponents;
    disjointSets.gatherComponents(2, rawComponents);

    // Store the components.
    components2.clear();
    for(const vector<uint64_t>& rawComponent: rawComponents) {
        Component& component2 = components2.emplace_back();
        component2.reserve(rawComponent.size());
        for(const uint64_t i: rawComponent) {
            component2.push_back(OrientedReadId::fromValue(ReadId(i)));
        }
    }
    components2.fillComponentId(orientedReadCount);

    cout << "Found " << components2.size() <<
        " non-trivial connected components of the read graph:" << endl;
    for(uint64_t componentId=0; componentId<components2.size(); componentId++) {
        const Component& component2 = components2[componentId];
        cout << "Component " << componentId << " has " <<
            component2.size() << " oriented reads and is " <<
            (component2.isDoubleStranded() ? "double" : "single") <<
            "-stranded." << endl;
    }
}



bool ReadGraph::Component::isDoubleStranded() const
{
    const Component& component = *this;

    if(size() < 2) {
        return false;
    } else {
        return component[0].getReadId() == component[1].getReadId();
    }
}



bool ReadGraph::Component::isSingleStranded() const
{
    return not isDoubleStranded();
}


void ReadGraph::Components::fillComponentId(uint64_t orientedReadCount)
{
    componentId.clear();
    componentId.resize(orientedReadCount, invalid<uint64_t>);
    for(uint64_t i=0; i<size(); i++) {
        const Component& component = (*this)[i];
        for(const OrientedReadId orientedReadId: component) {
            componentId[orientedReadId.getValue()] = i;
        }
    }
}



// We compute a single-stranded version of the ReadGraph as follows:
// - For each reverse complemented pair of single-stranded connected
//   component, we only keep one (the one in which the lowest numbered
//   OrientedReadId is on strand 0).
// - For each double-stranded component, we do a strand-aware
//   approximate min-cut to approximately separate strands.
//   The cut separates that component in two single-stranded
//   components and we keep only one of them (the one in which
//   the lowest numbered OrientedReadId is on strand 0).
void ReadGraph::computeComponents1()
{
    components1.clear();

    // Create a DisjointSets object that will be used repeatedly below.
    const ReadId readCount = anchors.reads.readCount();
    const uint64_t orientedReadCount = 2 * readCount;
    DisjointSets disjointSets(orientedReadCount);

    // Random generator used for shuffles.
    std::mt19937 random;

    // Loop over component of the complete double-stranded ReadGraph.
    for(uint64_t component2Id=0; component2Id<components2.size(); component2Id++) {
        const Component& component2 = components2[component2Id];
        const bool isDoubleStranded = component2.isDoubleStranded();
        cout << "Working on " << (isDoubleStranded ? "double" : "single") <<
            "-stranded component " << component2Id << " with " <<
            component2.size() << " orientedReads." << endl;

        // If a single stranded component, add it to components1
        // if its lowest number OrientedReaId is on strand0, and otherwise
        // do nothing with it.
        if(not isDoubleStranded) {
            if(component2.front().getStrand() == 0) {
                cout << "Copying this single-stranded component to the single-stranded ReadGraph." << endl;
            } else {
                cout << "Ignoring this single-stranded component." << endl;
            }
        } else {



            // This is a double-stranded component.
            // Compute an approximate, strand-aware min-cut using Karger's algorithm.
            cout << "Computing an approximate, strand-aware min-cut for this component." << endl;

            // Gather the indexes of EdgePairs in this component.
            vector<uint64_t> edgePairIndexes;
            for(uint64_t edgePairIndex=0; edgePairIndex<edgePairs.size(); edgePairIndex++) {
                const EdgePair& edgePair = edgePairs[edgePairIndex];
                if(edgePair.isCrossStrand) {
                    continue;
                }

                const OrientedReadId orientedReadId0(edgePair.readId0, 0);
                const OrientedReadId orientedReadId1(edgePair.readId1, edgePair.isSameStrand ? 0 : 1);
                uint64_t componentId = components2.componentId[orientedReadId0.getValue()];
                SHASTA2_ASSERT(componentId == components2.componentId[orientedReadId1.getValue()]);

                const OrientedReadId orientedReadId0rc(edgePair.readId0, 1);
                const OrientedReadId orientedReadId1rc(edgePair.readId1, edgePair.isSameStrand ? 1 : 0);
                uint64_t componentIdrc = components2.componentId[orientedReadId0rc.getValue()];
                SHASTA2_ASSERT(componentIdrc == components2.componentId[orientedReadId1rc.getValue()]);

                // Sanity check valid because we know we are in a double-stranded component.
                SHASTA2_ASSERT((componentId == component2Id) == (componentIdrc == component2Id));

                if(componentId == component2Id) {
                    edgePairIndexes.push_back(edgePairIndex);
                }
            }
            cout << "This component has " << 2 * edgePairIndexes.size() << " edges." << endl;



            // Iterate the Karger min-cut algorithm.
            vector<uint64_t> cutEdgePairIndexes;
            vector<uint64_t> bestCutEdgePairIndexes;
            vector< vector<uint64_t> > rawComponents;
            vector< vector<uint64_t> > bestRawComponents;
            uint64_t bestCutSize = std::numeric_limits<uint64_t>::max();
            for(uint64_t iteration=0; iteration<maxKargerIterationCount; iteration++) {
                cout << "Begin min-cut iteration " << iteration << endl;

                // Disconnect all OrientedReadIds in this component.
                for(const OrientedReadId orientedReadId: component2) {
                    disjointSets.initializeDisconnected(orientedReadId.getValue());
                }

                // Shuffle the edge pair indexes.
                std::shuffle(edgePairIndexes.begin(), edgePairIndexes.end(), random);

                // Process the EdgePairs in this shuffle order.
                cutEdgePairIndexes.clear();
                for(const uint64_t edgePairIndex: edgePairIndexes) {
                    const EdgePair& edgePair = edgePairs[edgePairIndex];

                    // Gather the OrientedReadIds involved in this EdgePair.
                    const OrientedReadId orientedReadId0(edgePair.readId0, 0);
                    const OrientedReadId orientedReadId1(edgePair.readId1, edgePair.isSameStrand ? 0 : 1);
                    const OrientedReadId orientedReadId0rc(edgePair.readId0, 1);
                    const OrientedReadId orientedReadId1rc(edgePair.readId1, edgePair.isSameStrand ? 1 : 0);

                    // Gather the corresponding components for the current cut.
                    const uint64_t component0 = disjointSets.findSet(orientedReadId0.getValue());
                    const uint64_t component1 = disjointSets.findSet(orientedReadId1.getValue());
                    const uint64_t component0rc = disjointSets.findSet(orientedReadId0rc.getValue());
                    const uint64_t component1rc = disjointSets.findSet(orientedReadId1rc.getValue());

                    // Sanity check that we are maintaining strand symmetry.
                    SHASTA2_ASSERT((component0 == component1) == (component0rc == component1rc));

                    // Sanity check that we are keeping strands separate.
                    SHASTA2_ASSERT(component0 != component0rc);
                    SHASTA2_ASSERT(component1 != component1rc);

                    // Only continue if not the same components
                    if(component0 != component1) {
                        SHASTA2_ASSERT(component0rc != component1rc);

                        // Find out if adding this adge pair would connect strands.
                        const bool wouldConnectStrands = (component0 == component1rc);
                        SHASTA2_ASSERT(wouldConnectStrands == (component0rc == component1));

                        if(wouldConnectStrands) {
                            cutEdgePairIndexes.push_back(edgePairIndex);
                        } else {
                            disjointSets.link(component0, component1);
                            disjointSets.link(component0rc, component1rc);
                        }
                    }
                }

                disjointSets.gatherComponents(2, rawComponents);

                cout << "Component " << component2Id << " iteration " << iteration <<
                    " produced a cut with " << 2 * cutEdgePairIndexes.size() << " edges and " <<
                    rawComponents.size() << " components." << endl;

                const uint64_t cutSize = cutEdgePairIndexes.size();
                if((rawComponents.size() == 2) and (cutSize < bestCutSize)) {
                    cout << "Updating the best cut." << endl;
                    bestCutSize = cutSize;
                    bestCutEdgePairIndexes = cutEdgePairIndexes;
                    bestRawComponents = rawComponents;
                }
            }

            cout << "The best cut for component " << component2Id <<
                " has " << 2 * bestCutSize << " edges." << endl;

            // Keep only one of the two components defined by the best cut.
            SHASTA2_ASSERT(bestRawComponents.size() == 2);
            const vector<uint64_t>& bestRawComponent0 = bestRawComponents[0];
            const vector<uint64_t>& bestRawComponent1 = bestRawComponents[1];
            const vector<uint64_t>& keepRawComponent =
                (bestRawComponent0.front() < bestRawComponent1.front()) ?
                bestRawComponent0 : bestRawComponent1;

            Component& component1 = components1.emplace_back();
            component1.reserve(keepRawComponent.size());
            for(const uint64_t i: keepRawComponent) {
                component1.push_back(OrientedReadId::fromValue(ReadId(i)));
            }
        }
    }

    components1.fillComponentId(orientedReadCount);
}



// Experiments with local strand-aware min-cuts.
void ReadGraph::localMinCut(
    ReadId readId,
    vector<uint64_t>& cutEdgePairIndexes) const
{
    // Gather OrientedReadIds that are at distance 1 from
    // readId-0 or readId-1.
    const span<const uint64_t> edgePairIndexesSpan = connectivityTable[readId];
    vector<uint64_t> edgePairIndexes(edgePairIndexesSpan.size());
    std::ranges::copy(edgePairIndexesSpan, edgePairIndexes.begin());

    vector<OrientedReadId> orientedReadIds;
    for(const uint64_t edgePairIndex: edgePairIndexes) {
        const EdgePair& edgePair = edgePairs[edgePairIndex];
        orientedReadIds.emplace_back(edgePair.readId0, 0);
        orientedReadIds.emplace_back(edgePair.readId0, 1);
        orientedReadIds.emplace_back(edgePair.readId1, 0);
        orientedReadIds.emplace_back(edgePair.readId1, 1);
    }
    deduplicate(orientedReadIds);
    cout << "The neighborhood of read " << readId <<
        " has " << orientedReadIds.size() << " vertices and " <<
        2 * edgePairIndexes.size() << " edges." << endl;

    cout << "The vertices correspond to the following OrientedReadIds:" << endl;
    for(uint64_t i=0; i<orientedReadIds.size(); i++) {
        cout << i << " " << orientedReadIds[i] << endl;
    }


    // Try a Karger mean cut.
    DisjointSets disjointSets(orientedReadIds.size());

    // Shuffle the edge pair indexes.
    std::mt19937 random;
    std::shuffle(edgePairIndexes.begin(), edgePairIndexes.end(), random);

    // Process the edge pairs in this order.
    cutEdgePairIndexes.clear();
    for(const uint64_t edgePairIndex: edgePairIndexes) {
        const EdgePair& edgePair = edgePairs[edgePairIndex];

        // Gather the OrientedReadIds involved in this EdgePair.
        const OrientedReadId orientedReadId0(edgePair.readId0, 0);
        const OrientedReadId orientedReadId1(edgePair.readId1, edgePair.isSameStrand ? 0 : 1);
        const OrientedReadId orientedReadId0rc(edgePair.readId0, 1);
        const OrientedReadId orientedReadId1rc(edgePair.readId1, edgePair.isSameStrand ? 1 : 0);

        cout << "Processing edge pair " << orientedReadId0 << " " << orientedReadId1 << " " <<
            orientedReadId0rc << " " << orientedReadId1rc << endl;

        // Gather the corresponding indexes in the OrientedReadIds.
        const auto it0 = std::ranges::lower_bound(orientedReadIds, orientedReadId0);
        const auto it1 = std::ranges::lower_bound(orientedReadIds, orientedReadId1);
        const auto it0rc = std::ranges::lower_bound(orientedReadIds, orientedReadId0rc);
        const auto it1rc = std::ranges::lower_bound(orientedReadIds, orientedReadId1rc);
        SHASTA2_ASSERT(it0 != orientedReadIds.end());
        SHASTA2_ASSERT(it1 != orientedReadIds.end());
        SHASTA2_ASSERT(it0rc != orientedReadIds.end());
        SHASTA2_ASSERT(it1rc != orientedReadIds.end());
        SHASTA2_ASSERT(*it0 == orientedReadId0);
        SHASTA2_ASSERT(*it1 == orientedReadId1);
        SHASTA2_ASSERT(*it0rc == orientedReadId0rc);
        SHASTA2_ASSERT(*it1rc == orientedReadId1rc);
        const uint64_t i0 = it0 - orientedReadIds.begin();
        const uint64_t i1 = it1 - orientedReadIds.begin();
        const uint64_t i0rc = it0rc - orientedReadIds.begin();
        const uint64_t i1rc = it1rc - orientedReadIds.begin();

        // Gather the corresponding components for the current cut.
        const uint64_t component0 = disjointSets.findSet(i0);
        const uint64_t component1 = disjointSets.findSet(i1);
        const uint64_t component0rc = disjointSets.findSet(i0rc);
        const uint64_t component1rc = disjointSets.findSet(i1rc);

        // Sanity check that we are maintaining strand symmetry.
        SHASTA2_ASSERT((component0 == component1) == (component0rc == component1rc));

        // Sanity check that we are keeping strands separate.
        SHASTA2_ASSERT(component0 != component0rc);
        SHASTA2_ASSERT(component1 != component1rc);

        // Only continue if not the same components
        if(component0 != component1) {
            SHASTA2_ASSERT(component0rc != component1rc);

            // Find out if adding this edge pair would connect strands.
            const bool wouldConnectStrands = (component0 == component1rc);
            SHASTA2_ASSERT(wouldConnectStrands == (component1 == component0rc));

            if(wouldConnectStrands) {
                cout << "Edge not added, becomes a cut edge." << endl;
                cutEdgePairIndexes.push_back(edgePairIndex);
            } else {
                cout << "Edge added." << endl;
                cout << "Merging " << component0 << " with " <<component1 << endl;
                cout << "Merging " << component0rc<< " with " <<component1rc << endl;
                disjointSets.link(component0, component1);
                disjointSets.link(component0rc, component1rc);
            }
        }


        // Big sanity check. REMOVE WHEN DONE DEBUGGING.
        SHASTA2_ASSERT((orientedReadIds.size() % 2) == 0);
        for(uint64_t i0=0; i0<orientedReadIds.size(); i0+=2) {
            const uint64_t i1 = i0 + 1;
            const OrientedReadId orientedReadId0 = orientedReadIds[i0];
            const OrientedReadId orientedReadId1 = orientedReadIds[i1];
            SHASTA2_ASSERT(orientedReadId0.getReadId() == orientedReadId1.getReadId());
            SHASTA2_ASSERT(orientedReadId0.getStrand() == 0);
            SHASTA2_ASSERT(orientedReadId1.getStrand() == 1);
            const uint64_t component0 = disjointSets.findSet(i0);
            const uint64_t component1 = disjointSets.findSet(i1);
            if(component0 == component1) {
                cout << "Assertion fails for " << orientedReadId0 << " " << orientedReadId1 << " " <<
                    component0 << " " << component1 << endl;
                SHASTA2_ASSERT(0);
            }
        }
    }

    cout << "Found a cut with " << cutEdgePairIndexes.size() << " edges." << endl;
    std::ranges::sort(cutEdgePairIndexes);

    vector< vector<uint64_t> > rawComponents;
    disjointSets.gatherComponents(1, rawComponents);
    cout << "This cut generates " << rawComponents.size() << " components with sizes";
    for(const vector<uint64_t>& rawComponent: rawComponents) {
        cout << " " << rawComponent.size();
    }
    cout << endl;


    // Graphviz output
    ofstream dot("LocalMinCut.dot");
    dot << "graph LocalMinCut {\n";

    for(const OrientedReadId orientedReadId: orientedReadIds) {
        const uint64_t i = std::ranges::lower_bound(orientedReadIds, orientedReadId) - orientedReadIds.begin();
        const bool isInComponent0 = std::ranges::binary_search(rawComponents[0], i);
        dot << "\"" << orientedReadId << "\" [color=" << (isInComponent0 == 0 ? "red" : "green") <<
            "];\n";
    }

    for(const uint64_t edgePairIndex: edgePairIndexes) {
        const EdgePair& edgePair = edgePairs[edgePairIndex];

        const bool isMinCutEdge = std::ranges::binary_search(cutEdgePairIndexes, edgePairIndex);

        const OrientedReadId orientedReadId0(edgePair.readId0, 0);
        const OrientedReadId orientedReadId1(edgePair.readId1, edgePair.isSameStrand ? 0 : 1);
        const OrientedReadId orientedReadId0rc(edgePair.readId0, 1);
        const OrientedReadId orientedReadId1rc(edgePair.readId1, edgePair.isSameStrand ? 1 : 0);

        dot << "\"" << orientedReadId0 << "\"--\"" << orientedReadId1 << "\"";
        if(isMinCutEdge) {
            dot << " [color=yellow]";
        }
        dot << ";\n";
        dot << "\"" << orientedReadId0rc << "\"--\"" << orientedReadId1rc << "\"";
        if(isMinCutEdge) {
            dot << " [color=yellow]";
        }
        dot << ";\n";
    }

    dot << "}\n";
}



void ReadGraph::writeConnectivityTable() const
{
    const ReadId readCount = anchors.reads.readCount();

    ofstream csv("ReadGraph-ConnectivityTable.csv");
    for(ReadId readId0=0; readId0<readCount; readId0++) {
        csv << readId0 << ",";

        const span<const uint64_t>& edgePairIndexes = connectivityTable[readId0];
        for(const uint64_t edgePairIndex: edgePairIndexes) {
            const EdgePair& edgePair = edgePairs[edgePairIndex];
            const ReadId readId1 = edgePair.getOther(readId0);
            csv << readId1 << " " << int(edgePair.isSameStrand) << ",";
        }

        csv << "\n";
    }

}



// Use self-complementary paths of length 2 to flag cross-strand EdgePairs.
void ReadGraph::flagCrossStrandEdgePairs2()
{
    const ReadId readCount = anchors.reads.readCount();
    vector<uint64_t> newCrossStrandEdgePairs;

    // Loop over ReadIds, looking for self-complementary quadrilaterals.
    for(ReadId readId=0; readId<readCount; readId++) {
        const span<const uint64_t>& edgePairIndexes = connectivityTable[readId];

        // Look for EdgePairs with the same other readId and opposite strands.
        for(uint64_t i1=1; i1<edgePairIndexes.size(); i1++) {
            const uint64_t i0 = i1 - 1;
            const uint64_t edgePairIndex0 = edgePairIndexes[i0];
            const uint64_t edgePairIndex1 = edgePairIndexes[i1];
            EdgePair& edgePair0 = edgePairs[edgePairIndex0];
            if(edgePair0.isCrossStrand) {
                continue;
            }
            EdgePair& edgePair1 = edgePairs[edgePairIndex1];
            if(edgePair1.isCrossStrand) {
                continue;
            }
            if(edgePair0.getOther(readId) == edgePair1.getOther(readId)) {
                SHASTA2_ASSERT(not edgePair0.isSameStrand);
                SHASTA2_ASSERT(edgePair1.isSameStrand);
                newCrossStrandEdgePairs.push_back(edgePairIndex0);
                newCrossStrandEdgePairs.push_back(edgePairIndex1);
            }
        }
    }


    ofstream csv("EdgePairs2.csv");
    deduplicate(newCrossStrandEdgePairs);
    for(uint64_t edgePairIndex: newCrossStrandEdgePairs) {
        EdgePair& edgePair = edgePairs[edgePairIndex];
        csv << edgePair.readId0 << "," << edgePair.readId1 << endl;
    }
    cout <<
        "Of " << 2 * edgePairs.size() << " read graph edges, " <<
        2 * newCrossStrandEdgePairs.size() <<
        " are in self-complementary paths of length 2." << endl;

}



// This can be multithreaded.
void ReadGraph::findStrandSymmetriQuadrilaterals()
{
    const ReadId readCount = anchors.reads.readCount();
    strandSymmetricQuadrilaterals.clear();

    // Loop over ReadIds.
    for(ReadId readId0=0; readId0<readCount; readId0++) {
        const span<const uint64_t>& edgePairIndexes = connectivityTable[readId0];

        // Loop over EdgePairs that this ReadId is involved in,
        // looking for EdgePairs with the same other readId and opposite strands.
        for(uint64_t i1=1; i1<edgePairIndexes.size(); i1++) {
            const uint64_t i0 = i1 - 1;

            const uint64_t edgePairIndex0 = edgePairIndexes[i0];
            const uint64_t edgePairIndex1 = edgePairIndexes[i1];

            EdgePair& edgePair0 = edgePairs[edgePairIndex0];
            SHASTA2_ASSERT(edgePair0.readId0 < edgePair0.readId1);

            EdgePair& edgePair1 = edgePairs[edgePairIndex1];
            SHASTA2_ASSERT(edgePair1.readId0 < edgePair1.readId1);


            if(edgePair0.getOther(readId0) == edgePair1.getOther(readId0)) {
                SHASTA2_ASSERT(not edgePair0.isSameStrand);
                SHASTA2_ASSERT(edgePair1.isSameStrand);
                const ReadId readId1 = edgePair0.getOther(readId0);

                // These two EdgePairs defined a StrandSymmetricQuadrilateral.
                // It will be found twice, so only store it if readId0 < readId1.
                if(readId0 < readId1) {
                    StrandSymmetricQuadrilateral& strandSymmetricQuadrilateral =
                        strandSymmetricQuadrilaterals.emplace_back();
                    strandSymmetricQuadrilateral.readId0 = readId0;
                    strandSymmetricQuadrilateral.readId1 = readId1;
                    strandSymmetricQuadrilateral.edgePairIndexSameStrand = edgePairIndex1;
                    strandSymmetricQuadrilateral.edgePairIndexOppositeStrands = edgePairIndex0;
                    strandSymmetricQuadrilateral.coverageSameStrand = edgePair1.coverage;
                    strandSymmetricQuadrilateral.coverageOppositeStrands = edgePair0.coverage;
                }
            }
        }
    }

    cout << "Found " << strandSymmetricQuadrilaterals.size() <<
        " strand-symmetric quadrilaterals in the ReadGraph." << endl;

    {
        ofstream csv("StrandSymmetricQuadrilaterals.csv");
        csv << "ReadId0,ReadId1,CoverageSameStrand,CoverageOppositeStrands,\n";
        for(const auto& strandSymmetricQuadrilateral: strandSymmetricQuadrilaterals) {
            csv << strandSymmetricQuadrilateral.readId0 << ",";
            csv << strandSymmetricQuadrilateral.readId1 << ",";
            csv << strandSymmetricQuadrilateral.coverageSameStrand << ",";
            csv << strandSymmetricQuadrilateral.coverageOppositeStrands << ",";
            csv << "\n";;
        }
    }


    // Count the number of quadrilaterals each readId is involved in.
    vector<uint64_t> histogram(readCount, 0);
    for(const auto& strandSymmetricQuadrilateral: strandSymmetricQuadrilaterals) {
        ++histogram[strandSymmetricQuadrilateral.readId0];
        ++histogram[strandSymmetricQuadrilateral.readId1];
    }

    {
        ofstream csv("StrandSymmetricQuadrilaterals-Histogram.csv");
        csv << "ReadId,Count,\n";
        for(ReadId readId=0; readId<readCount; readId++) {
            csv << readId << ",";
            csv << histogram[readId] << ",";
            csv << "\n";
        }
    }
}



// Use self-complementary paths of length m to flag cross-strand EdgePairs.
uint64_t ReadGraph::flagCrossStrandEdgePairs(uint64_t maxDistance)
{
    const ReadId readCount = anchors.reads.readCount();

    // vector<uint64_t> newCrossStrandEdgePairs;
    uint64_t newCrossStrandEdgePairsCount = 0;

    // Loop over ReadIds.
    for(ReadId readId=0; readId<readCount; readId++) {
        const OrientedReadId orientedReadIdA(readId, 0);
        const OrientedReadId orientedReadIdB(readId, 1);

        // Look for a path of length noo more than m
        // between orientedReadId0 and orientedReadId1.
        // This is a quick and dirty version.
        // A more efficient implementation is possible.

        std::queue<OrientedReadId> q;
        q.push(orientedReadIdA);

        // Information about the vertices we already encountered.
        class VertexInfo {
        public:
            uint64_t parentEdgePairIndex = invalid<uint64_t>;
            uint64_t distance = invalid<uint64_t>;
            VertexInfo() {}
            VertexInfo(uint64_t parentEdgePairIndex, uint64_t distance) :
                parentEdgePairIndex(parentEdgePairIndex),
                distance(distance) {}
        };
        std::map<OrientedReadId, VertexInfo> vertexMap;
        vertexMap.insert(make_pair(orientedReadIdA, VertexInfo(invalid<uint64_t>, 0)));



        // BFS loop.
        // cout << "BFS loop beginning at " <<
        bool found = false;
        while(not q.empty()) {
            const OrientedReadId orientedReadId0 = q.front();
            q.pop();
            const ReadId readId0 = orientedReadId0.getReadId();
            const Strand strand0 = orientedReadId0.getStrand();
            const uint64_t distance0 = vertexMap[orientedReadId0].distance;
            const uint64_t distance1 = distance0 + 1;

            // Loop over neighbors of orientedReadId0.
            const span<const uint64_t> edgePairIndexes = connectivityTable[readId0];
            for(uint64_t edgePairIndex: edgePairIndexes) {
                const EdgePair& edgePair = edgePairs[edgePairIndex];
                if(edgePair.isCrossStrand) {
                    continue;
                }
                const ReadId readId1 = edgePair.getOther(readId0);
                const OrientedReadId orientedReadId1(readId1, edgePair.isSameStrand ? strand0 : 1 - strand0);

                const auto it1 = vertexMap.find(orientedReadId1);
                if(it1 != vertexMap.end()) {
                    continue;
                }
                vertexMap.insert(make_pair(orientedReadId1, VertexInfo(edgePairIndex, distance1)));
                if(orientedReadId1 == orientedReadIdB) {
                    found = true;
                    break;
                }

                if(distance1 < maxDistance) {
                    q.push(orientedReadId1);
                }
            }
        }

        if(found) {

            if(true) {
                cout << "Found a path of length " << maxDistance << " for " << readId <<endl;
                cout << "Reverse path is:";
                OrientedReadId orientedReadId = orientedReadIdB;
                while(true) {
                    cout << " " << orientedReadId;
                    if(orientedReadId == orientedReadIdA) {
                        break;
                    }
                    const auto it = vertexMap.find(orientedReadId);
                    SHASTA2_ASSERT(it != vertexMap.end());
                    const uint64_t edgePairIndex = it->second.parentEdgePairIndex;
                    const EdgePair& edgePair = edgePairs[edgePairIndex];
                    cout << " " << edgePair.coverage;
                    const ReadId otherReadId = edgePair.getOther(orientedReadId.getReadId());
                    const Strand otherStrand = (edgePair.isSameStrand ? orientedReadId.getStrand() : 1 - orientedReadId.getStrand());
                    orientedReadId = OrientedReadId(otherReadId, otherStrand);
                }
                cout << endl;
            }


#if 0
            // Add to newCrossStrandEdgePairs all the edge pair indexes on this path.
            OrientedReadId orientedReadId = orientedReadIdB;
            while(true) {
                if(orientedReadId == orientedReadIdA) {
                    break;
                }
                const auto it = vertexMap.find(orientedReadId);
                SHASTA2_ASSERT(it != vertexMap.end());
                const uint64_t edgePairIndex = it->second.parentEdgePairIndex;
                newCrossStrandEdgePairs.push_back(edgePairIndex);
                const EdgePair& edgePair = edgePairs[edgePairIndex];
                const ReadId otherReadId = edgePair.getOther(orientedReadId.getReadId());
                const Strand otherStrand = (edgePair.isSameStrand ? orientedReadId.getStrand() : 1 - orientedReadId.getStrand());
                orientedReadId = OrientedReadId(otherReadId, otherStrand);
            }
#endif

#if 0

            // Add to newCrossStrandEdgePairs the EdgePair with lowest coverage on this path.
            uint64_t lowestCoverageEdgePairIndex = invalid<uint64_t>;
            uint64_t lowestCoverage = std::numeric_limits<uint64_t>::max();
            OrientedReadId orientedReadId = orientedReadIdB;
            while(true) {
                if(orientedReadId == orientedReadIdA) {
                    break;
                }
                const auto it = vertexMap.find(orientedReadId);
                SHASTA2_ASSERT(it != vertexMap.end());
                const uint64_t edgePairIndex = it->second.parentEdgePairIndex;
                const EdgePair& edgePair = edgePairs[edgePairIndex];
                if(edgePair.coverage < lowestCoverage) {
                    lowestCoverage = edgePair.coverage;
                    lowestCoverageEdgePairIndex = edgePairIndex;
                }
                const ReadId otherReadId = edgePair.getOther(orientedReadId.getReadId());
                const Strand otherStrand = (edgePair.isSameStrand ? orientedReadId.getStrand() : 1 - orientedReadId.getStrand());
                orientedReadId = OrientedReadId(otherReadId, otherStrand);
            }
            // newCrossStrandEdgePairs.push_back(lowestCoverageEdgePairIndex);
            edgePairs[lowestCoverageEdgePairIndex].isCrossStrand = true;
            ++newCrossStrandEdgePairsCount;
#endif
        }
    }

#if 0
    deduplicate(newCrossStrandEdgePairs);
    for(uint64_t edgePairIndex: newCrossStrandEdgePairs) {
        EdgePair& edgePair = edgePairs[edgePairIndex];
        edgePair.isCrossStrand = true;
    }
#endif
    cout <<
        "Of " << 2 * edgePairs.size() << " read graph edges, " <<
        2 * newCrossStrandEdgePairsCount <<
        " where flagged as cross-strand using self-complementary paths of length " << maxDistance << endl;

    return newCrossStrandEdgePairsCount;
}

