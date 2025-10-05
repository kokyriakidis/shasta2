#include "AssemblyGraph.hpp"
#include "Anchor.hpp"
#include "Assembler.hpp"
#include "deduplicate.hpp"
#include "RestrictedAnchorGraph.hpp"
#include "TangleMatrix1.hpp"
using namespace shasta;


// Constructor from another AssemblyGraph and assembly paths
// like the ones computed by ReadFollowing.
// (Those are not actually paths in the AssemblyGraph but simply
// sequences of edges).
AssemblyGraph::AssemblyGraph(
    const Assembler& assembler,
    const Options& options,
    const AssemblyGraph& oldAssemblyGraph,
    const vector< vector<edge_descriptor> >& assemblyPaths) :
    MappedMemoryOwner(assembler),
    MultithreadedObject<AssemblyGraph>(*this),
    anchors(assembler.anchors()),
    journeys(assembler.journeys()),
    options(options),
    orderById(*this)
{
    AssemblyGraph& newAssemblyGraph = *this;

    // The terminal edges of the paths are special in that
    // only one copy of them is carried in the new AssemblyGraph.

    // Gather the terminal edges.
    vector<edge_descriptor> terminalEdges;
    for(const vector<edge_descriptor>& path: assemblyPaths) {
        terminalEdges.push_back(path.front());
        terminalEdges.push_back(path.back());
    }
    deduplicate(terminalEdges);
    sort(terminalEdges.begin(), terminalEdges.end(), oldAssemblyGraph.orderById);



    // Generate an edge in the newAssemblyGraph for each terminal
    // edge in the old assembly graph.
    // Keep track of them, keyed by edge_descriptor in the oldAssemblyGraph.
    std::map<edge_descriptor, edge_descriptor> terminalEdgesMap;
    for(const edge_descriptor eOld: terminalEdges) {
        const AssemblyGraphEdge& oldEdge = oldAssemblyGraph[eOld];
        const AnchorId anchorId0 = oldEdge.front().anchorPair.anchorIdA;
        const AnchorId anchorId1 = oldEdge.back().anchorPair.anchorIdB;
        const vertex_descriptor v0 = add_vertex(AssemblyGraphVertex(anchorId0, newAssemblyGraph.nextVertexId++), newAssemblyGraph);
        const vertex_descriptor v1 = add_vertex(AssemblyGraphVertex(anchorId1, newAssemblyGraph.nextVertexId++), newAssemblyGraph);
        edge_descriptor eNew;
        tie(eNew, ignore) = add_edge(v0, v1, AssemblyGraphEdge(oldEdge.id), newAssemblyGraph);
        AssemblyGraphEdge& newEdge = newAssemblyGraph[eNew];
        newEdge = oldEdge;
        newEdge.id = newAssemblyGraph.nextEdgeId++;
        terminalEdgesMap.insert(make_pair(eOld, eNew));
    }



    // Now generate the remaining edges for each path.
    // Note that non-terminal edges can get more than one copy
    // in the newAssemblyGraph.
    vector< vector<edge_descriptor> > newAssemblyPaths;
    for(const auto& oldAssemblyPath: assemblyPaths) {
        newAssemblyPaths.emplace_back();
        auto& newAssemblyPath = newAssemblyPaths.back();

        for(const edge_descriptor eOld: oldAssemblyPath) {
            const auto it = terminalEdgesMap.find(eOld);
            if(it != terminalEdgesMap.end()) {
                newAssemblyPath.push_back(it->second);
            } else {

                // This is not a terminal edge and we have to create.
                const AssemblyGraphEdge& oldEdge = oldAssemblyGraph[eOld];
                const AnchorId anchorId0 = oldEdge.front().anchorPair.anchorIdA;
                const AnchorId anchorId1 = oldEdge.back().anchorPair.anchorIdB;
                const vertex_descriptor v0 = add_vertex(AssemblyGraphVertex(anchorId0, newAssemblyGraph.nextVertexId++), newAssemblyGraph);
                const vertex_descriptor v1 = add_vertex(AssemblyGraphVertex(anchorId1, newAssemblyGraph.nextVertexId++), newAssemblyGraph);
                edge_descriptor eNew;
                tie(eNew, ignore) = add_edge(v0, v1, AssemblyGraphEdge(oldEdge.id), newAssemblyGraph);
                AssemblyGraphEdge& newEdge = newAssemblyGraph[eNew];
                newEdge = oldEdge;
                newEdge.id = newAssemblyGraph.nextEdgeId++;
                newAssemblyPath.push_back(eNew);
            }
        }
    }



    // For each pair of consecutive edges in a path (e0, e1),
    // generate a new edge in-between to bridge between them.
    // The code is similar to Tangle1::addConnectPair and Tangle1::detangle,
    // but simpler.
    for(const auto& newAssemblyPath: newAssemblyPaths) {

        for(uint64_t i1=1; i1<newAssemblyPath.size(); i1++) {
            const uint64_t i0 = i1 - 1;
            const edge_descriptor e0 = newAssemblyPath[i0];
            const edge_descriptor e1 = newAssemblyPath[i1];

            const vertex_descriptor v0 = target(e0, newAssemblyGraph);
            const vertex_descriptor v1 = source(e1, newAssemblyGraph);

            const AnchorId anchorId0 = newAssemblyGraph[v0].anchorId;
            const AnchorId anchorId1 = newAssemblyGraph[v1].anchorId;

            // Create the new edge.
            // If the two anchors are the same, leave it empty without any steps.
            // Otherwise use the same process in Tangle1::addConnectPair.
            edge_descriptor eNew;
            tie(eNew, ignore) = add_edge(v0, v1, AssemblyGraphEdge(newAssemblyGraph.nextEdgeId++), newAssemblyGraph);
            AssemblyGraphEdge& newEdge = newAssemblyGraph[eNew];
            if(anchorId0 != anchorId1) {

                // Create the RestrictedAnchorGraph, then:
                // - Remove vertices not accessible from anchorId0 and anchorId1.
                // - Remove cycles.
                // - Find the longest path.
                // - Add one step for each edge of the longest path of the RestrictedAnchorGraph.

                ostream html(0);
                const TangleMatrix1 tangleMatrix(
                    newAssemblyGraph,
                    vector<edge_descriptor>(1, e0),
                    vector<edge_descriptor>(1, e1),
                    html);

                RestrictedAnchorGraph restrictedAnchorGraph(
                    newAssemblyGraph.anchors, newAssemblyGraph.journeys, tangleMatrix, 0, 0, html);
                restrictedAnchorGraph.removeLowCoverageEdges(anchorId0, anchorId1);
                restrictedAnchorGraph.keepBetween(anchorId0, anchorId1);
                restrictedAnchorGraph.removeCycles();
                restrictedAnchorGraph.keepBetween(anchorId0, anchorId1);
                vector<RestrictedAnchorGraph::edge_descriptor> longestPath;
                // restrictedAnchorGraph.findLongestPath(longestPath);
                restrictedAnchorGraph.findOptimalPath(anchorId0, anchorId1, longestPath);

                for(const RestrictedAnchorGraph::edge_descriptor re: longestPath) {
                    const auto& rEdge = restrictedAnchorGraph[re];
                    if(rEdge.anchorPair.size() < newAssemblyGraph.options.detangleMinCoverage) {
                        newEdge.clear();
                        SHASTA_ASSERT(0);
                    }
                    newEdge.push_back(AssemblyGraphEdgeStep(rEdge.anchorPair, rEdge.offset));
                }
            }
        }
    }

    newAssemblyGraph.removeEmptyEdges();
    newAssemblyGraph.compress();
}
