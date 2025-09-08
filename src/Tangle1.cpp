// Shasta.
#include "Tangle1.hpp"
#include "Anchor.hpp"
#include "RestrictedAnchorGraph.hpp"
#include "TangleMatrix1.hpp"
using namespace shasta;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>



// Constructor from a set of AssemblyGraph vertices.
Tangle1::Tangle1(
    AssemblyGraph& assemblyGraph,
    const vector<vertex_descriptor>& tangleVerticesArgument) :
    assemblyGraph(assemblyGraph),
    tangleVertices(tangleVerticesArgument)
{
    // Sort the tangle vertices by id.
    std::ranges::sort(tangleVertices, assemblyGraph.orderById);

    // Find the edges.
    findEntrances();
    findExits();
    findTangleEdges();

    // Compute the tangle matrix.
    ostream html(0);
    tangleMatrixPointer = make_shared<TangleMatrix1>(assemblyGraph, entrances, exits, html);
}



// Constructor from a single vertex.
Tangle1::Tangle1(
    AssemblyGraph& assemblyGraph,
    vertex_descriptor v) :
    Tangle1(
        assemblyGraph,
        vector<vertex_descriptor>(1, v))
{}



// Constructor from an edge.
Tangle1::Tangle1(
    AssemblyGraph& assemblyGraph,
    edge_descriptor e) :
    Tangle1(
        assemblyGraph,
        vector<vertex_descriptor>({source(e, assemblyGraph), target(e, assemblyGraph)}))
{}



void Tangle1::findEntrances()
{
    entrances.clear();
    for(const vertex_descriptor v1: tangleVertices) {
        BGL_FORALL_INEDGES(v1, e, assemblyGraph, AssemblyGraph) {
            const vertex_descriptor v0 = source(e, assemblyGraph);
            if(not isTangleVertex(v0)) {
                entrances.push_back(e);
            }
        }
    }
    std::ranges::sort(entrances, assemblyGraph.orderById);
}



void Tangle1::findExits()
{
    exits.clear();
    for(const vertex_descriptor v0: tangleVertices) {
        BGL_FORALL_OUTEDGES(v0, e, assemblyGraph, AssemblyGraph) {
            const vertex_descriptor v1 = target(e, assemblyGraph);
            if(not isTangleVertex(v1)) {
                exits.push_back(e);
            }
        }
    }
    std::ranges::sort(exits, assemblyGraph.orderById);
}



void Tangle1::findTangleEdges()
{
    tangleEdges.clear();
    for(const vertex_descriptor v0: tangleVertices) {
        BGL_FORALL_OUTEDGES(v0, e, assemblyGraph, AssemblyGraph) {
            const vertex_descriptor v1 = target(e, assemblyGraph);
            if(isTangleVertex(v1)) {
                tangleEdges.push_back(e);
            }
        }
    }
}



bool Tangle1::isTangleVertex(vertex_descriptor v) const
{
    return std::ranges::binary_search(tangleVertices, v, assemblyGraph.orderById);
}



void Tangle1::addConnectPair(uint64_t entranceIndex, uint64_t exitIndex) {
    SHASTA_ASSERT(entranceIndex < entrances.size());
    SHASTA_ASSERT(exitIndex < exits.size());
    connectPairs.emplace_back(entranceIndex, exitIndex);
    AssemblyGraphEdge& newEdge = connectPairs.back().newEdge;

    // Store the AssemblyGraphEdge, without adding it to the AssemblyGraph for now.
    const AssemblyGraph::vertex_descriptor v0 = target(entrances[entranceIndex], assemblyGraph);
    const AssemblyGraph::vertex_descriptor v1 = source(exits[exitIndex], assemblyGraph);
    const AnchorId anchorId0 = assemblyGraph[v0].anchorId;
    const AnchorId anchorId1 = assemblyGraph[v1].anchorId;

    if(anchorId0 == anchorId1) {
        // We just generate an empty edge without any steps.
        // If any of these are left after compress, they will have to be
        // removed by collapsing the vertices they join.
    } else {
        // Create the RestrictedAnchorGraph, then:
        // - Remove vertices not accessible from anchorId0 and anchorId1.
        // - Remove cycles.
        // - Find the longest path.
        // - Add one step for each edge of the longest path of the RestrictedAnchorGraph.
        ostream html(0);
        RestrictedAnchorGraph restrictedAnchorGraph(
            assemblyGraph.anchors, assemblyGraph.journeys, tangleMatrix(), entranceIndex, exitIndex, html);
        restrictedAnchorGraph.removeLowCoverageEdges(anchorId0, anchorId1);
        restrictedAnchorGraph.keepBetween(anchorId0, anchorId1);
        restrictedAnchorGraph.removeCycles();
        restrictedAnchorGraph.keepBetween(anchorId0, anchorId1);
        vector<RestrictedAnchorGraph::edge_descriptor> longestPath;
        // restrictedAnchorGraph.findLongestPath(longestPath);
        restrictedAnchorGraph.findOptimalPath(anchorId0, anchorId1, longestPath);

        for(const RestrictedAnchorGraph::edge_descriptor re: longestPath) {
            const auto& rEdge = restrictedAnchorGraph[re];
            newEdge.push_back(AssemblyGraphEdgeStep(rEdge.anchorPair,rEdge.offset));
        }
    }
}



void Tangle1::detangle()
{
    // Reroute the entrances to new vertices, so all
    // entrances and exits become temporarily dangling.
    vector<vertex_descriptor> newEntranceVertices;
    rerouteEntrances(newEntranceVertices);

    vector<vertex_descriptor> newExitVertices;
    rerouteExits(newExitVertices);

    // Finally, reconnect the entrance/exit pairs in our ConnectPairs.
    for(ConnectPair& connectPair: connectPairs) {
        const uint64_t entranceIndex = connectPair.entranceIndex;
        const uint64_t exitIndex = connectPair.exitIndex;
        reconnect(
            connectPair,
            newEntranceVertices[entranceIndex],
            newExitVertices[exitIndex]);
    }

    // Now we can remove all the tangle vertices and their edges.
    // This also removes the old entrances and exits.
    removedVertices = tangleVertices;
    for(const vertex_descriptor v: tangleVertices) {
        if(false) {
            cout << "Removing for detangling:";
            BGL_FORALL_INEDGES(v, e, assemblyGraph, AssemblyGraph) {
                cout << " " << assemblyGraph[e].id;
            }
            BGL_FORALL_OUTEDGES(v, e, assemblyGraph, AssemblyGraph) {
                cout << " " << assemblyGraph[e].id;
            }
            cout << endl;
        }
        clear_vertex(v, assemblyGraph);
        remove_vertex(v, assemblyGraph);
    }
}



// Make a copy of each entrance edge, with the target vertex replaced by a new vertex
// with the same AnchorId.
void Tangle1::rerouteEntrances(vector<vertex_descriptor>& newEntranceVertices) const
{
    newEntranceVertices.clear();

    for(const edge_descriptor& eOld: entrances) {
        const vertex_descriptor v0Old = source(eOld, assemblyGraph);
        AssemblyGraphEdge& edgeOld = assemblyGraph[eOld];
        const AnchorId lastAnchorId = edgeOld.back().anchorPair.anchorIdB;

        // Create the new vertex.
        const vertex_descriptor v1 = add_vertex(
            AssemblyGraphVertex(lastAnchorId, assemblyGraph.nextVertexId++), assemblyGraph);
        newEntranceVertices.push_back(v1);

        // Create the new edge, with the same steps and id as the old one.
        edge_descriptor eNew;
        tie(eNew, ignore) = add_edge(v0Old, v1, AssemblyGraphEdge(edgeOld.id), assemblyGraph);
        AssemblyGraphEdge& edgeNew = assemblyGraph[eNew];
        edgeNew.swapSteps(edgeOld);
        edgeNew.wasAssembled = edgeOld.wasAssembled;
    }
}



// Make a copy of each exit edge, with the source vertex replaced by a new vertex
// with the same AnchorId.
void Tangle1::rerouteExits(vector<vertex_descriptor>& newExitVertices) const
{
    newExitVertices.clear();
    for(const edge_descriptor& eOld: exits) {
        const vertex_descriptor v1Old = target(eOld, assemblyGraph);
        AssemblyGraphEdge& edgeOld = assemblyGraph[eOld];
        const AnchorId firstAnchorId = edgeOld.front().anchorPair.anchorIdA;

        // Create the new vertex.
        const vertex_descriptor v0 = add_vertex(
            AssemblyGraphVertex(firstAnchorId, assemblyGraph.nextVertexId++), assemblyGraph);
        newExitVertices.push_back(v0);

        // Create the new edge, with the same steps and id as the old one.
        edge_descriptor eNew;
        tie(eNew, ignore) = add_edge(v0, v1Old, AssemblyGraphEdge(edgeOld.id), assemblyGraph);
        AssemblyGraphEdge& edgeNew = assemblyGraph[eNew];
        edgeNew.swapSteps(edgeOld);
        edgeNew.wasAssembled = edgeOld.wasAssembled;
    }
}



void Tangle1::reconnect(
    ConnectPair& connectPair,
    vertex_descriptor v0,
    vertex_descriptor v1
    ) const
{
    const bool debug = false;

    const AnchorId anchorId0 = assemblyGraph[v0].anchorId;
    const AnchorId anchorId1 = assemblyGraph[v1].anchorId;

    if(debug) {
        cout << "Connecting entrance " << assemblyGraph[entrances[connectPair.entranceIndex]].id <<
            " with exit " << assemblyGraph[exits[connectPair.exitIndex]].id << endl;
        cout << "Anchors " << anchorIdToString(anchorId0) << " " <<
            anchorIdToString(anchorId1) << endl;
    }


    // Create a new AssemblyGraphEdge between v0 and v1.
    AssemblyGraph::edge_descriptor e;
    tie(e, ignore) = add_edge(v0, v1, assemblyGraph);
    AssemblyGraphEdge& newEdge = assemblyGraph[e];
    newEdge.swapSteps(connectPair.newEdge);
    newEdge.id = assemblyGraph.nextEdgeId++;
    if(debug) {
        cout << "Created new assembly graph edge " << newEdge.id << endl;
    }


#if 0
    if(anchorId0 == anchorId1) {
        // We just generate an empty edge without any steps.
        // If any of these are left after compress, they will have to be
        // removed by collapsing the vertices they join.
    } else {
        // Create the RestrictedAnchorGraph, then:
        // - Remove vertices not accessible from anchorId0 and anchorId1.
        // - Remove cycles.
        // - Find the longest path.
        // - Add one step for each edge of the longest path of the RestrictedAnchorGraph.
        ostream html(0);
        RestrictedAnchorGraph restrictedAnchorGraph(
            assemblyGraph.anchors, assemblyGraph.journeys, tangleMatrix(), iEntrance, iExit, html);
        restrictedAnchorGraph.keepBetween(anchorId0, anchorId1);
        restrictedAnchorGraph.removeCycles();
        restrictedAnchorGraph.keepBetween(anchorId0, anchorId1);
        vector<RestrictedAnchorGraph::edge_descriptor> longestPath;
        // restrictedAnchorGraph.findLongestPath(longestPath);
        restrictedAnchorGraph.findOptimalPath(anchorId0, anchorId1, longestPath);

        for(const RestrictedAnchorGraph::edge_descriptor re: longestPath) {
            const auto& rEdge = restrictedAnchorGraph[re];
            newEdge.push_back(AssemblyGraphEdgeStep(rEdge.anchorPair,rEdge.offset));
        }
    }
#endif
}



uint64_t Tangle1::minConnectPairCoverage() const
{
    uint64_t minCoverage = std::numeric_limits<uint64_t>::max();
    for(const ConnectPair& connectPair: connectPairs) {
        for(const auto& step: connectPair.newEdge) {
            minCoverage = min(minCoverage, step.anchorPair.size());
        }
    }
    return minCoverage;
}
