// Shasta.
#include "StrandSeparator1.hpp"
#include "MarkerKmers.hpp"
#include "Reads.hpp"
#include "timestamp.hpp"
using namespace shasta2;


StrandSeparator1::StrandSeparator1(
    const Anchors& anchors2,
    Anchors& anchors1
    ) :
    anchors2(anchors2),
    anchorCount(anchors2.anchorMarkerInfos.size()),
    readCount(anchors2.reads.readCount()),
    orientedReadCount(2 * readCount),
    vertexCount(anchorCount + orientedReadCount),
    rank(vertexCount),
    parent(vertexCount),
    disjointSets(&rank[0], &parent[0])
{
    cout << timestamp << "Strand separation begins." << endl;
    cout << "Number of reads is " << readCount << endl;
    cout << "Number of oriented reads is " << orientedReadCount << endl;
    cout << "Total number of anchors before strand separation is " << anchorCount << endl;
    SHASTA2_ASSERT((anchorCount % 2) == 0);

    gatherEdgePairs();
    computeComponents2();
    createComponents1();
    generateSingleStrandedAnchors(anchors1);

    cout << timestamp << "Strand separation ends." << endl;

}



void StrandSeparator1::gatherEdgePairs()
{
    edgePairs.clear();

    // Loop over positive anchors (the AnchorId is even).
    for(AnchorId anchorId=0; anchorId<anchorCount; anchorId+=2) {
        const Anchor anchor = anchors2[anchorId];
        for(const AnchorMarkerInfo& anchorMarkerInfo: anchor) {
            const OrientedReadId orientedReadId = anchorMarkerInfo.orientedReadId;
            const ReadId readId = orientedReadId.getReadId();
            const Strand strand = orientedReadId.getStrand();
            const bool onSameStrand = (strand == 0);

            edgePairs.push_back(EdgePair(anchorId, readId, onSameStrand));
        }
    }

    cout << "The double-stranded bipartite graph has " <<
        edgePairs.size() << " edge pairs." << endl;

}



void StrandSeparator1::initializeDisjointSets()
{
    for(uint64_t vertexId=0; vertexId<vertexCount; vertexId++) {
        disjointSets.make_set(vertexId);
    }
}



void StrandSeparator1::computeComponents2()
{
    initializeDisjointSets();
    for(const EdgePair& edgePair: edgePairs) {
        applyEdgePairToDisjointSets(edgePair);
    }
    createComponents(components2);

    cout << "Found " << components2.components.size() <<
        " connected components of the double-stranded bipartite graph." << endl;
}



void StrandSeparator1::applyEdgePairToDisjointSets(const EdgePair& edgePair)
{
    const VertexIds vertexIds =  getVertexIds(edgePair);
    disjointSets.union_set(vertexIds.anchorVertexA, vertexIds.orientedReadVertexA);
    disjointSets.union_set(vertexIds.anchorVertexB, vertexIds.orientedReadVertexB);
}



void StrandSeparator1::Component::check() const
{
    SHASTA2_ASSERT(not isTrivial());
    SHASTA2_ASSERT(std::ranges::is_sorted(anchorIds));
    SHASTA2_ASSERT(std::ranges::is_sorted(orientedReadIds));
}



void StrandSeparator1::Components::initialize(uint64_t anchorCount, uint64_t orientedReadCount)
{
    components.clear();
    anchorComponent.clear();
    orientedReadComponent.clear();
    anchorComponent.resize(anchorCount, invalid<uint64_t>);
    orientedReadComponent.resize(orientedReadCount, invalid<uint64_t>);
}



void StrandSeparator1::Components::fill()
{
    for(uint64_t componentId=0; componentId<components.size(); componentId++) {
        const Component& component = components[componentId];
        component.check();

        // Fill anchorComponent.
        for(const AnchorId anchorId: component.anchorIds) {
            uint64_t& thisAnchorComponent = anchorComponent[anchorId];
            SHASTA2_ASSERT(thisAnchorComponent == invalid<uint64_t>);
            thisAnchorComponent = componentId;
        }

        // Fill orientedReadComponent.
        for(OrientedReadId orientedReadId: component.orientedReadIds) {
            uint64_t& thisOrientedReadComponent = orientedReadComponent[orientedReadId.getValue()];
            SHASTA2_ASSERT(thisOrientedReadComponent == invalid<uint64_t>);
            thisOrientedReadComponent = componentId;
        }
    }
}



void StrandSeparator1::createComponents(Components& components)
{

    vector< vector<uint64_t> > componentTable(vertexCount);
    for(uint64_t vertexId=0; vertexId<vertexCount; vertexId++) {
        const uint64_t componentId = disjointSets.find_set(vertexId);
        componentTable[componentId].push_back(vertexId);
    }

    components.initialize(anchorCount, orientedReadCount);
    for(const vector<uint64_t>& componentVertices: componentTable) {
        if(componentVertices.size() > 1) {
            components.components.emplace_back();
            Component& component = components.components.back();
            for(const uint64_t vertexId: componentVertices) {
                if(vertexId < anchorCount) {
                    const AnchorId anchorId = vertexId;
                    component.anchorIds.push_back(anchorId);
                } else {
                    const OrientedReadId orientedReadId =
                        OrientedReadId::fromValue(ReadId(vertexId - anchorCount));
                    component.orientedReadIds.push_back(orientedReadId);
                }
            }
        }
    }

    components.fill();
}



void StrandSeparator1::createComponents1()
{
    components1.initialize(anchorCount, orientedReadCount);

    initializeDisjointSets();
    for(uint64_t componentId2=0; componentId2<components2.components.size(); componentId2++) {
        createComponent1(componentId2);
    }

    components1.fill();
}



// This assumes that the disjointSets are initially all disjoint.
// It also restores the disjointSets to the completely disjoint state.
void StrandSeparator1::createComponent1(uint64_t componentId2)
{
    const Component& component2 = components2.components[componentId2];
    cout << "Working on component " << componentId2 <<
        " of the double-stranded bipartite graph." << endl;
    cout << "This component has " << component2.orientedReadIds.size() <<
        " oriented reads and " << component2.anchorIds.size() <<
        " anchors." << endl;
    component2.check();



    // If the component2 is single-stranded, it is one of a pair
    // of reverse complemented pair of single-stranded components.
    // We want to keep only one the two in the pair.
    // We keep the one in which the lowest number ReadId
    // is on strand 0.
    if(not component2.isDoubleStranded()) {
        cout << "This component is single-stranded." << endl;
        if(component2.orientedReadIds[0].getStrand() == 0) {
            components1.components.emplace_back(component2);
            cout << "This single-stranded component was stored." << endl;
        } else {
            cout << "This single-stranded component was not stored." << endl;
        }
        return;
    }


    // If getting here, this component2 is double-stranded and we have to
    // separate it.
    SHASTA2_ASSERT(component2.isDoubleStranded());
    cout << "This component is double-stranded and requires strand separation." << endl;



    // Gather the edge pair ids for edges in this connected component.
    vector<uint64_t> edgePairIds;
    // Loop over positive anchors in this component2.
    for(const AnchorId anchorId: component2.anchorIds) {
        if((anchorId % 2) == 1) {
            continue;
        }
        const Anchor anchor =  anchors2[anchorId];

        // Loop over AnchorMarkerInfos for this anchor.
        // Each of them generates an EdgePair.
        const uint64_t firstMarkerInfoId = &anchor.front() - anchors2.anchorMarkerInfos.begin();
        const uint64_t firstEdgePairId = firstMarkerInfoId / 2;
        for(uint64_t i=0; i<anchor.size(); i++) {
            const uint64_t edgePairId = firstEdgePairId + i;
            const EdgePair& edgePair = edgePairs[edgePairId];

            // Sanity checks on the EdgePair.
            SHASTA2_ASSERT(edgePair.anchorId == anchorId);
            SHASTA2_ASSERT(edgePair.anchorIdA() == anchorId);
            SHASTA2_ASSERT(edgePair.anchorIdB() == (anchorId ^ 1UL));
            SHASTA2_ASSERT(components2.anchorComponent[edgePair.anchorIdA()] == componentId2);
            SHASTA2_ASSERT(components2.anchorComponent[edgePair.anchorIdB()] == componentId2);
            SHASTA2_ASSERT(components2.orientedReadComponent[edgePair.orientedReadIdA().getValue()] == componentId2);
            SHASTA2_ASSERT(components2.orientedReadComponent[edgePair.orientedReadIdB().getValue()] == componentId2);

            // Store this EdgePair.
            edgePairIds.push_back(edgePairId);
        }
    }
    cout << "In the bipartite graph, this double-stranded component has " <<
        edgePairIds.size() << " edge pairs." << endl;



    // Iterate. At each iteration we generate a cut between strands
    // using a strand-aware version of Krager's min-cut algorithm.
    const uint64_t iterationCount = 5000;
    vector<uint64_t> bestCutEdgePairIds;
    bool hadSuccessfulIteration = false;
    Component bestComponent1;
    for(uint64_t iteration=0; iteration<iterationCount; iteration++) {
        cout << timestamp << "Starting strand separation iteration " << iteration <<
            " for this connected component." << endl;

        // Remove all edges from the counterpart of this component in the
        // single-stranded bipartite graph.
        for(const AnchorId anchorId: component2.anchorIds) {
            const uint64_t anchorVertex = anchorId;
            disjointSets.make_set(anchorVertex);
        }
        for(const OrientedReadId orientedReadId: component2.orientedReadIds) {
            const uint64_t orientedReadVertex = anchorCount + orientedReadId.getValue();
            disjointSets.make_set(orientedReadVertex);
        }

        // Shuffle the edge pair ids.
        std::shuffle(edgePairIds.begin(), edgePairIds.end(), random);



        // Add the edge pairs in shuffled order to the single-stranded bipartite graph,
        // but don't add the ones that would cause strand mixing.
        vector<uint64_t> cutEdgePairIds;
        for(const uint64_t edgePairId: edgePairIds) {
            const EdgePair& edgePair = edgePairs[edgePairId];

            // Get the VertexIds for this EdgePair.
            const VertexIds vertexIds = getVertexIds(edgePair);

            // Get the corresponding components in the single-stranded bipartite graph.
            const uint64_t anchorComponentA = disjointSets.find_set(vertexIds.anchorVertexA);
            const uint64_t anchorComponentB = disjointSets.find_set(vertexIds.anchorVertexB);
            const uint64_t orientedReadComponentA = disjointSets.find_set(vertexIds.orientedReadVertexA);
            const uint64_t orientedReadComponentB = disjointSets.find_set(vertexIds.orientedReadVertexB);

            // Check if adding this edge would mix strand.
            if(anchorComponentA == orientedReadComponentB) {
                SHASTA2_ASSERT(anchorComponentB == orientedReadComponentA);
                // Adding this edge would cause strand mixing. Don't add the edge.
                cutEdgePairIds.push_back(edgePairId);
            } else {
                SHASTA2_ASSERT(anchorComponentB != orientedReadComponentA);
                if(anchorComponentA != orientedReadComponentA) {
                    disjointSets.link(anchorComponentA, orientedReadComponentA);
                }
                if(anchorComponentB != orientedReadComponentB) {
                    disjointSets.link(anchorComponentB, orientedReadComponentB);
                }
            }
        }
        cout << "The cut at this iteration contains " << cutEdgePairIds.size() <<
            " edge pairs." << endl;

        // Figure out how many single-stranded components we got.
        Components components1;
        createComponents(components1);
        cout << "The cut at this iteration generated " << components1.components.size() <<
            " single-stranded components." << endl;

        // Check that they are all indeed single-stranded.
        for(const Component& component1: components1.components) {
            SHASTA2_ASSERT(not component1.isDoubleStranded());
        }


#if 0
        // To find strand contact regions, compute connected components using
        // only the edges of the best cut.
        // Remove all edges from the single-stranded counterpart of this component.
        cout << "Looking for strand contact regions given this cut." << endl;
        for(const AnchorId anchorId: component2.anchorIds) {
            const uint64_t anchorVertex = anchorId;
            disjointSets.make_set(anchorVertex);
        }
        for(const OrientedReadId orientedReadId: component2.orientedReadIds) {
            const uint64_t orientedReadVertex = anchorCount + orientedReadId.getValue();
            disjointSets.make_set(orientedReadVertex);
        }
        for(const uint64_t edgePairId: cutEdgePairIds) {
            const EdgePair& edgePair = edgePairs[edgePairId];
            applyEdgePairToDisjointSets(edgePair);
        }
        Components contactRegions;
        createComponents(contactRegions);
        cout << "Found " << contactRegions.components.size() <<
            " strand contact regions:" << endl;
        for(const Component& component: contactRegions.components) {
            cout << component.anchorIds.size() << " anchors, " <<
                component.orientedReadIds.size() << " oriented reads, " <<
                (component.isDoubleStranded() ? "double" : "single") << "-stranded." << endl;
        }
#endif


        // Only use this cut if it generated exactly 2 single-stranded connected components.
        if(components1.components.size() == 2) {
            if((not hadSuccessfulIteration) or (cutEdgePairIds.size() < bestCutEdgePairIds.size())) {
                hadSuccessfulIteration = true;
                cout << "Updating the best cut." << endl;
                bestCutEdgePairIds = cutEdgePairIds;
                bestComponent1 = components1.components.front();
            }
        }
    }

    // Store the best single-stranded component we found.
    components1.components.emplace_back(bestComponent1);
}



void StrandSeparator1::generateSingleStrandedAnchors(Anchors& anchors1) const
{
    // Initialize the kmerToAnchorTable.
    anchors1.kmerToAnchorTable.resize(anchors2.markerKmers.size());
    std::ranges::fill(anchors1.kmerToAnchorTable, invalid<AnchorId>);


    // Loop over connected components of the single-stranded bipartite graph.
    for(uint64_t componentId1=0; componentId1<components1.components.size(); componentId1++) {
        const Component& component1 = components1.components[componentId1];

        // Loop over anchors in this component.
        // Each Anchor in this component generates a new single-stranded anchor.
        for(const AnchorId anchorId2: component1.anchorIds) {
            const Anchor anchor2 = anchors2[anchorId2];

            // Get the anchorId for the new anchor we are generating.
            const AnchorId anchorId1 = anchors1.size();

            // Copy the AnchorInfo.
            const AnchorInfo& anchorInfo2 = anchors2.anchorInfos[anchorId2];
            anchors1.anchorInfos.push_back(anchorInfo2);

            // Update the kmerToAnchorTable.
            anchors1.kmerToAnchorTable[anchorInfo2.kmerIndex] = anchorId1;

            // Generate the AnchorMarkerInfos for this new anchor.
            // They are the same as for anchor2, but we can only use
            // OrientedReadIds in the same component.
            anchors1.anchorMarkerInfos.appendVector();
            for(const AnchorMarkerInfo& anchorMarkerInfo: anchor2) {
                const OrientedReadId orientedReadId= anchorMarkerInfo.orientedReadId;
                if(components1.orientedReadComponent[orientedReadId.getValue()] == componentId1) {
                    anchors1.anchorMarkerInfos.append(anchorMarkerInfo);
                }
            }

        }
    }

}
