// Shasta.
#include "AssemblyGraph.hpp"
#include "Anchor.hpp"
#include "AnchorGraph.hpp"
#include "deduplicate.hpp"
#include "findLinearChains.hpp"
using namespace shasta;

// Standard library.
#include <fstream.hpp>
#include <queue>
#include <tuple.hpp>



AssemblyGraph::AssemblyGraph(
    const Anchors& anchors,
    const AnchorGraph& anchorGraph)
{
    AssemblyGraph& assemblyGraph = *this;

    // Find linear chains in the AnchorGraph.
    vector< vector<AnchorGraph::edge_descriptor> > chains;
    findLinearChains(anchorGraph, 0, chains);
    cout << "Found " << chains.size() << " linear chains in the anchor graph." << endl;

    // Get the initial and final AnchorId of all chains.
    // These generate an AssemblyGraph vertex each, after deduplication.
    vector<AnchorId> vertexAnchorIds;
    for(const auto& chain: chains) {
        const AnchorGraph::edge_descriptor eA = chain.front();
        const AnchorGraph::edge_descriptor eB = chain.back();
        const AnchorId anchorIdA = anchorGraph[eA].anchorIdA;
        const AnchorId anchorIdB = anchorGraph[eB].anchorIdB;
        vertexAnchorIds.push_back(anchorIdA);
        vertexAnchorIds.push_back(anchorIdB);
    }
    deduplicate(vertexAnchorIds);

    // Create the vertices.
    std::map<AnchorId, vertex_descriptor> vertexMap;
    for(const AnchorId anchorId: vertexAnchorIds) {
        const vertex_descriptor v = add_vertex(AssemblyGraphVertex(anchorId), assemblyGraph);
        vertexMap.insert(make_pair(anchorId, v));
    }



    // Create the edges.
    for(const auto& chain: chains) {
        const AnchorGraph::edge_descriptor eA = chain.front();
        const AnchorGraph::edge_descriptor eB = chain.back();
        const AnchorId anchorIdA = anchorGraph[eA].anchorIdA;
        const AnchorId anchorIdB = anchorGraph[eB].anchorIdB;
        const vertex_descriptor vA = vertexMap[anchorIdA];
        const vertex_descriptor vB = vertexMap[anchorIdB];

        edge_descriptor e;
        tie(e, ignore) = add_edge(vA, vB, assemblyGraph);
        AssemblyGraphEdge& edge = assemblyGraph[e];
        edge.id = nextEdgeId++;

        for(const AnchorGraph::edge_descriptor e: chain) {
            edge.resize(edge.size() + 1);
            AssemblyGraphStep& assemblyGraphStep = edge.back();
            assemblyGraphStep.anchorPair = anchorGraph[e];
            assemblyGraphStep.anchorPair.getOffsetStatistics(
                anchors,
                assemblyGraphStep.averageOffset,
                assemblyGraphStep.minOffset,
                assemblyGraphStep.maxOffset);
        }
        edge.computeOffsets();
    }

    cout << "The initial assembly graph has " << num_vertices(assemblyGraph) <<
        " vertices and " << num_edges(assemblyGraph) << " edges." << endl;
    write("A");

    transitiveReduction();
    compress();

    cout << "After transitive reduction, the assembly graph has " << num_vertices(assemblyGraph) <<
        " vertices and " << num_edges(assemblyGraph) << " edges." << endl;
    write("B");
}



void AssemblyGraphEdge::computeOffsets()
{
    averageOffset = 0;
    minOffset = 0;
    maxOffset = 0;
    for(const AssemblyGraphStep& step: *this) {
        averageOffset += step.averageOffset;
        minOffset += step.minOffset;
        maxOffset += step.maxOffset;
    }
}



void AssemblyGraph::writeGfa(const string& fileName) const
{
    const AssemblyGraph& assemblyGraph = *this;

    ofstream gfa(fileName);

    // Write the header line.
    gfa << "H\tVN:Z:1.0\n";

    // Write a segment for each edge.
    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        const AssemblyGraphEdge& edge = assemblyGraph[e];

        // Record type.
        gfa << "S\t";

        // Name.
        gfa << edge.id << "\t";

        // Sequence.
        gfa << "*\t";

        // Sequence length in bases.
        gfa << "LN:i:" << edge.averageOffset << "\n";
    }

    // For each vertex, write links between each pair of incoming/outgoing edges.
    BGL_FORALL_VERTICES(v, assemblyGraph, AssemblyGraph) {
        BGL_FORALL_INEDGES(v, e0, assemblyGraph, AssemblyGraph) {
            const uint64_t id0 = assemblyGraph[e0].id;
            BGL_FORALL_OUTEDGES(v, e1, assemblyGraph, AssemblyGraph) {
                const uint64_t id1 = assemblyGraph[e1].id;
                gfa <<
                    "L\t" <<
                    id0 << "\t+\t" <<
                    id1 << "\t+\t*\n";
            }
        }
    }
}



void AssemblyGraph::write(const string& name) const
{
    writeGfa("AssemblyGraph-" + name + ".gfa");
    writeSegments("AssemblyGraphSegments-" + name + ".csv");
    writeSegmentDetails("AssemblyGraphSegmentsDetails-" + name + ".csv");
}



void AssemblyGraph::writeSegments(const string& fileName) const
{
    const AssemblyGraph& assemblyGraph = *this;

    ofstream csv(fileName);
    csv << "Id,AnchorIdA,AnchorIdB,StepCount,AverageOffset,MinOffset,MaxOffset\n";

    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        const AssemblyGraphEdge& edge = assemblyGraph[e];

        csv <<
            edge.id << "," <<
            anchorIdToString(edge.anchorIdA()) << "," <<
            anchorIdToString(edge.anchorIdB()) << "," <<
            edge.size() << "," <<
            edge.averageOffset << "," <<
            edge.minOffset << "," <<
            edge.maxOffset << "\n";
    }

}



void AssemblyGraph::writeSegmentDetails(const string& fileName) const
{
    const AssemblyGraph& assemblyGraph = *this;

    ofstream csv(fileName);
    csv << "Id,Position,AnchorIdA,AnchorIdB,AverageOffset,MinOffset,MaxOffset\n";

    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        const AssemblyGraphEdge& edge = assemblyGraph[e];

        for(uint64_t i=0; i<edge.size(); i++) {
            const AssemblyGraphStep& step = edge[i];

            csv <<
                edge.id << "," <<
                i << "," <<
                anchorIdToString(step.anchorPair.anchorIdA) << "," <<
                anchorIdToString(step.anchorPair.anchorIdB) << "," <<
                step.averageOffset << "," <<
                step.minOffset << "," <<
                step.maxOffset << "\n";
        }
    }

}



// An v0->v1 edge with length length01 is removed if:
// - length01 is no more than threshold1.
// - A BFS that starts at v0 and does not use edge v0->v1
//   reaches v1 with a BFS offset no more than a  + b * length01.
// Most of the edges removed in this way are caused by errors.
// Even if a correct one is removed at this edge, it can still
// be assembled correctly by the LocalAssembly, is phasing is
// successful in this region.
void AssemblyGraph::transitiveReduction()
{
    AssemblyGraph& assemblyGraph = *this;

    // EXPOSE WHEN CODE STABILIZES.
    const uint64_t threshold1 = 10000;
    const uint64_t a = 200;
    const uint64_t b = 2;

    // Work areas for the BFS are defined here to reduce nemory allocation activity.
    std::queue<vertex_descriptor> q;
    vector<vertex_descriptor> verticesEncountered;

    // Loop over all edges. We will be removing edges while iterating,
    // so we have to avoid invalidating iterators.
    edge_iterator it, itEnd;
    tie(it, itEnd) = edges(assemblyGraph);
    while(it != itEnd) {
        const edge_descriptor e = *it;
        const AssemblyGraphEdge& edge = assemblyGraph[e];

        // Increment the iterator here, before possibly removing this edge.
        ++it;

        if(edge.length() > threshold1) {
            continue;
        }

        // This edge is a candidate for removal.
        // Initialize the BFS.
        const vertex_descriptor v0 = source(e, assemblyGraph);
        const vertex_descriptor v1 = target(e, assemblyGraph);
        SHASTA_ASSERT(q.empty());
        SHASTA_ASSERT(verticesEncountered.empty());
        q.push(v0);
        verticesEncountered.push_back(v0);
        assemblyGraph[v0].bfsDistance = 0;
        const uint64_t threshold2 = a + b * edge.length();

        // Main BFS loop.
        bool removeEdge = false;
        while(not q.empty()) {
            const vertex_descriptor vA = q.front();
            q.pop();
            const uint64_t distanceA = assemblyGraph[vA].bfsDistance;

            // Loop over out-edges of vA.
            BGL_FORALL_OUTEDGES(vA, eAB, assemblyGraph, AssemblyGraph) {

                // Don't use edge e in the BFS.
                if(eAB == e) {
                    continue;
                }

                const vertex_descriptor vB = target(eAB, assemblyGraph);

                // If we already encountered vB, do nothing.
                if(assemblyGraph[vB].bfsDistance != invalid<uint64_t>) {
                    continue;
                }

                // If we got too far, don't use this in the BFS.
                const uint64_t distanceB = distanceA + assemblyGraph[eAB].length();
                if(distanceB > threshold2) {
                    continue;
                }

                // If vB is v1, edge e should be removed.
                if(vB == v1) {
                    removeEdge = true;
                    break;
                }

                // Update the BFS.
                q.push(vB);
                verticesEncountered.push_back(vB);
                assemblyGraph[vB].bfsDistance = distanceB;

            }

            if(removeEdge) {
                break;
            }
        }

        if(removeEdge) {
            boost::remove_edge(e, assemblyGraph);
        }

        // Cleanup the BFS.
        for(const vertex_descriptor v: verticesEncountered) {
            assemblyGraph[v].bfsDistance = invalid<uint64_t>;
        }
        while(not q.empty()) {
            q.pop();
        }
        verticesEncountered.clear();
    }
}



void AssemblyGraph::compress()
{
    AssemblyGraph& assemblyGraph = *this;

    // Find linear chains. All the edges in a linear chain of
    // more than one edge are combined into a single edge.
    vector< vector<edge_descriptor> > linearChains;
    findLinearChains(assemblyGraph, 0, linearChains);
    cout << "Compress found " << linearChains.size() << " linear chains." << endl;

    for(const vector<edge_descriptor>& linearChain: linearChains) {
        if(linearChain.size() == 1) {
            continue;
        }

        const vertex_descriptor v0 = source(linearChain.front(), assemblyGraph);
        const vertex_descriptor v1 = target(linearChain.back(), assemblyGraph);

        // Create a new edge consisting of just this linear chain.
        edge_descriptor eNew;
        tie(eNew, ignore) = add_edge(v0, v1, assemblyGraph);
        AssemblyGraphEdge& newEdge = assemblyGraph[eNew];
        newEdge.id = nextEdgeId++;

        // Concatenate the old edges to create the new one,
        // then remove the old edges.
        newEdge.averageOffset = 0;
        newEdge.minOffset = 0;
        newEdge.maxOffset = 0;
        for(const edge_descriptor eOld: linearChain) {
            AssemblyGraphEdge& oldEdge = assemblyGraph[eOld];
            newEdge.averageOffset += oldEdge.averageOffset;
            newEdge.minOffset += oldEdge.minOffset;
            newEdge.maxOffset += oldEdge.maxOffset;
            copy(oldEdge.begin(), oldEdge.end(), back_inserter(newEdge));
            boost::remove_edge(eOld, assemblyGraph);
        }

        // Now we can remove the middle vertices in the linear chain.
        for(uint64_t i=1; i<linearChain.size(); i++) {
            const vertex_descriptor v = source(linearChain[i], assemblyGraph);
            remove_vertex(v, assemblyGraph);
        }
    }
}
