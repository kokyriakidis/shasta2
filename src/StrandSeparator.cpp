// Shasta.
#include "StrandSeparator.hpp"
#include "orderPairs.hpp"
#include "Reads.hpp"
#include "timestamp.hpp"
using namespace shasta2;

// Boost libraries.
#include <boost/pending/disjoint_sets.hpp>

// Standard library.
#include <random>


StrandSeparator::StrandSeparator(

    // The input double-stranded Anchors.
    // Anchors must be in consecutively numbered pairs of
    // reverse complemented anchors.
    const Anchors& anchors2,

    // The output single-stranded Anchors.
    // They must be empty on input.
    Anchors& /* anchors1 */
    ) :
    anchorCount(anchors2.anchorMarkerInfos.size()),
    readCount(anchors2.reads.readCount()),
    orientedReadCount(2 * readCount)
{
    cout << timestamp << "Strand separation begins." << endl;
    cout << "Number of reads is " << readCount << endl;
    cout << "Number of oriented reads is " << orientedReadCount << endl;
    cout << "Total number of anchors before strand separation is " << anchorCount << endl;

    gatherEdges(anchors2);
    computeComponents2();
    computeComponents1(anchors2);

    cout << timestamp << "Strand separation ends." << endl;
}



// Each AnchorMarkerInfo in anchors2 generates an Edge.
void StrandSeparator::gatherEdges(const Anchors& anchors2)
{
    edges.clear();

    for(AnchorId anchorId=0; anchorId<anchorCount; anchorId++) {
        const Anchor anchor = anchors2[anchorId];
        for(const AnchorMarkerInfo& anchorMarkerInfo: anchor) {
            edges.push_back(Edge(anchorId, anchorMarkerInfo.orientedReadId));
        }
    }

    cout << "The double-stranded bipartite graph has " <<
        edges.size() << " edges." << endl;
}



void StrandSeparator::computeComponents2()
{

    // We use a disjoint set with vertices numbered as follows:
    // - The vertex corresponding to AnchorId anchorId
    //   has index anchorId.
    // - The vertex corresponding to OrientedReadId orientedReadId
    //   has index anchorCount + orientedReadId.getValue().

    // Initialize the disjoint sets data structure.
    const uint64_t n = anchorCount + orientedReadCount;
    vector<uint64_t> rank(n);
    vector<uint64_t> parent(n);
    boost::disjoint_sets<uint64_t*, uint64_t*> disjointSets(&rank[0], &parent[0]);
    for(uint64_t i=0; i<n; i++) {
        disjointSets.make_set(i);
    }

    // Do a union_set operation for each Edge.
    for(const Edge& edge: edges) {
        const uint64_t anchorVertex = edge.anchorId;
        const uint64_t orientedReadVertex = anchorCount + edge.orientedReadId.getValue();
        disjointSets.union_set(anchorVertex, orientedReadVertex);
    }

    // Gather the connected components, as numbered by the disjointSets.
    vector<Component> rawComponents(n);
    for(AnchorId anchorId=0; anchorId<anchorCount; anchorId++) {
        const uint64_t anchorVertex = anchorId;
        const uint64_t rawComponentId = disjointSets.find_set(anchorVertex);
        rawComponents[rawComponentId].anchorIds.push_back(anchorId);
    }
    for(ReadId readId=0; readId<readCount; readId++) {
        for(Strand strand=0; strand<2; strand++) {
            const OrientedReadId orientedReadId(readId, strand);
            const uint64_t orientedReadVertex = anchorCount + orientedReadId.getValue();
            const uint64_t rawComponentId = disjointSets.find_set(orientedReadVertex);
            rawComponents[rawComponentId].orientedReadIds.push_back(orientedReadId);
        }
    }

    // Gather the rawComponentId of components with more than one oriented read
    // and store them sorted by decreasing number of oriented reads.
    vector< pair<uint64_t, uint64_t> > rawComponentTable;   // (rawComponentId, number of oriented reads)
    for(uint64_t rawComponentId=0; rawComponentId<n; rawComponentId++) {
        const Component& rawComponent = rawComponents[rawComponentId];
        if(rawComponent.orientedReadIds.size() > 1) {
            rawComponentTable.emplace_back(rawComponentId, rawComponent.orientedReadIds.size());
        }
    }
    sort(rawComponentTable.begin(), rawComponentTable.end(),
        OrderPairsBySecondOnlyGreater<uint64_t, uint64_t>());

    // Store the components with this numbering.
    components2.clear();
    for(const auto& p: rawComponentTable) {
        const uint64_t rawComponentId = p.first;
        const uint64_t orientedReadCount = p.second;
        const Component& rawComponent = rawComponents[rawComponentId];
        SHASTA2_ASSERT(rawComponent.orientedReadIds.size() == orientedReadCount);
        components2.emplace_back(rawComponent);
    }
    cout << "The double-stranded bipartite graph has " << components1.size() <<
        " connected components." << endl;

    // Store anchorComponent2 and orientedReadComponent2.
    anchorComponent2.clear();
    anchorComponent2.resize(anchorCount, invalid<uint64_t>);
    orientedReadComponent2.clear();
    orientedReadComponent2.resize(orientedReadCount, invalid<uint64_t>);
    for(uint64_t componentId2=0; componentId2<components2.size(); componentId2++) {
        const Component& component2 = components2[componentId2];
        for(const AnchorId anchorId: component2.anchorIds) {
            anchorComponent2[anchorId] = componentId2;
        }
        for(const OrientedReadId orientedReadId: component2.orientedReadIds) {
            orientedReadComponent2[orientedReadId.getValue()] = componentId2;
        }
    }
}



void StrandSeparator::computeComponents1(const Anchors& anchors2)
{
    // Work vector used below and defined here to reduce
    // memory allocation activity.
    vector<uint64_t> edgeIds;

    // Random generator used for shuffles.
    std::mt19937 random;

    // Initialize the disjoint sets data structure for the single-stranded biopartite graph.
    const uint64_t n = anchorCount + orientedReadCount;
    vector<uint64_t> rank(n);
    vector<uint64_t> parent(n);
    boost::disjoint_sets<uint64_t*, uint64_t*> disjointSets(&rank[0], &parent[0]);
    for(uint64_t i=0; i<n; i++) {
        disjointSets.make_set(i);
    }

    // Initialize the components of the single-stranded bipartite graph.
    components1.clear();



    // Process one component of the double-stranded bipartite graph at a time.
    for(uint64_t componentId2=0; componentId2<components2.size(); componentId2++) {
        const Component& component2 = components2[componentId2];
        cout << "Working on connected component " << componentId2 << " with " <<
        component2.orientedReadIds.size() << " oriented reads, " <<
            component2.anchorIds.size() << " anchors." << endl;
        SHASTA2_ASSERT(component2.orientedReadIds.size() >= 2);

        // Figure out if this is a double-stranded component
        const bool isDoubleStranded =
            component2.orientedReadIds[0].getReadId() == component2.orientedReadIds[1].getReadId();

        // Single-stranded component.
        // We want to keep only one of each pair of complementari single-stranded components.
        // So we choose the one where the lowest number OrientedReadId is on strand 0.
        if(not isDoubleStranded) {
            cout << "This connected component is single-stranded." << endl;
            if(component2.orientedReadIds.front().getStrand() == 0) {
                components1.emplace_back(component2);
            }
        }

        // Double-stranded component.
        else {
            cout << "This connected component is double-stranded." << endl;

            // Gather the edge ids for edges in this connected component.
            edgeIds.clear();
            for(const AnchorId anchorId: component2.anchorIds) {
                const Anchor anchor =  anchors2[anchorId];
                const uint64_t firstEdgeId = &anchor.front() - anchors2.anchorMarkerInfos.begin();
                for(uint64_t i=0; i<anchor.size(); i++) {
                    const AnchorMarkerInfo& anchorMarkerInfo = anchor[i];
                    SHASTA2_ASSERT(
                        orientedReadComponent2[anchorMarkerInfo.orientedReadId.getValue()] ==
                        componentId2);
                    edgeIds.push_back(firstEdgeId + i);
                }
            }
            cout << "In the bipartite graph, this double-stranded component has " <<
                edgeIds.size() << " edges." << endl;



            // Iterate. At each iteration we use a strand-aware version of
            // Krager's min-cut algorithm.
            const uint64_t iterationCount = 6;
            uint64_t bestSkippedCount = std::numeric_limits<uint64_t>::max();
            Component component1;
            for(uint64_t iteration=0; iteration<iterationCount; iteration++) {

                // Shuffle the edge ids.
                std::shuffle(edgeIds.begin(), edgeIds.end(), random);

                // Remove all edges from the single-stranded counterpart of this component.
                for(const AnchorId anchorId: component2.anchorIds) {
                    const uint64_t anchorVertex = anchorId;
                    disjointSets.make_set(anchorVertex);
                }
                for(const OrientedReadId orientedReadId: component2.orientedReadIds) {
                    const uint64_t orientedReadVertex = anchorCount + orientedReadId.getValue();
                    disjointSets.make_set(orientedReadVertex);
                }

                // Add the edges in shuffled order to the single-stranded bipartite graph,
                // but don't add the ones that would cause strand mixing.
                // Every time we add an edge we also add its reverse complement.
                uint64_t skippedCount = 0;
                for(const uint64_t edgeId: edgeIds) {
                    const Edge& edge = edges[edgeId];
                    const AnchorId anchorIdA = edge.anchorId;
                    const OrientedReadId orientedReadIdA = edge.orientedReadId;
                    OrientedReadId orientedReadIdB = orientedReadIdA;
                    orientedReadIdB.flipStrand();
                    const AnchorId anchorIdB = anchorIdA ^ 1UL;

                    // Get the corresponding vertices.
                    const uint64_t anchorVertexA = anchorIdA;
                    const uint64_t anchorVertexB = anchorIdB;
                    const uint64_t orientedReadVertexA = anchorCount + orientedReadIdA.getValue();
                    const uint64_t orientedReadVertexB = anchorCount + orientedReadIdB.getValue();

                    // Get the corresponding components.
                    const uint64_t anchorComponentA = disjointSets.find_set(anchorVertexA);
                    const uint64_t anchorComponentB = disjointSets.find_set(anchorVertexB);
                    const uint64_t orientedReadComponentA = disjointSets.find_set(orientedReadVertexA);
                    const uint64_t orientedReadComponentB = disjointSets.find_set(orientedReadVertexB);

                    // Check if adding this edge would mix strand.
                    if(anchorComponentA == orientedReadComponentB) {
                        SHASTA2_ASSERT(anchorComponentB == orientedReadComponentA);
                        // Adding this edge would cause strand mixing. Don't add the edge.
                        ++skippedCount;
                    } else {
                        SHASTA2_ASSERT(anchorComponentB != orientedReadComponentA);
                        disjointSets.union_set(anchorVertexA, orientedReadVertexA);
                        disjointSets.union_set(anchorVertexB, orientedReadVertexB);
                    }
                }
                cout << "Iteration " << iteration << ": skipped " << skippedCount << " edges to avoid strand mixing." << endl;

                // Figure out how many single-stranded components we ended up with.
                std::map<uint64_t, Component> component1Map;
                for(const AnchorId anchorId: component2.anchorIds) {
                    const uint64_t anchorVertex = anchorId;
                    const uint64_t anchorComponent = disjointSets.find_set(anchorVertex);
                    component1Map[anchorComponent].anchorIds.push_back(anchorId);
                }
                for(const OrientedReadId orientedReadId: component2.orientedReadIds) {
                    const uint64_t orientedReadVertex = anchorCount + orientedReadId.getValue();
                    const uint64_t orientedReadComponent = disjointSets.find_set(orientedReadVertex);
                    component1Map[orientedReadComponent].orientedReadIds.push_back(orientedReadId);
                }
                cout << "Found " << component1Map.size() << " single-stranded components." << endl;
                /*
                for(const auto& p: component1Map) {
                    const Component& component = p.second;
                    cout << component.orientedReadIds.size() << " oriented reads, " <<
                        component.anchorIds.size() << " anchors." << endl;
                }
                */

                // Only use it if we have exactly 2 single-stranded connected components.
                if(component1Map.size() == 2) {
                    if(skippedCount < bestSkippedCount) {
                        cout << "Updating the single-stranded component." << endl;
                        bestSkippedCount = skippedCount;
                        component1 = component1Map.begin()->second;
                    }
                } else {
                    cout << "This iteration did not create exactly 2 components." << endl;
                }
            }
            components1.emplace_back(component1);
            cout << "The single-stranded component has " <<
                component1.orientedReadIds.size() << " oriented reads and " <<
                component1.anchorIds.size() << " anchors." << endl;
        }
    }
    cout << "The single-stranded bipartite graph has " << components1.size() <<
        " connected components." << endl;



    // Store anchorComponent1 and orientedReadComponent1.
    anchorComponent1.clear();
    anchorComponent1.resize(anchorCount, invalid<uint64_t>);
    orientedReadComponent1.clear();
    orientedReadComponent1.resize(orientedReadCount, invalid<uint64_t>);
    for(uint64_t componentId1=0; componentId1<components1.size(); componentId1++) {
        const Component& component1 = components1[componentId1];
        for(const AnchorId anchorId: component1.anchorIds) {
            anchorComponent1[anchorId] = componentId1;
        }
        for(const OrientedReadId orientedReadId: component1.orientedReadIds) {
            orientedReadComponent1[orientedReadId.getValue()] = componentId1;
        }
    }
}
