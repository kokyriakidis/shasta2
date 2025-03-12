// Shasta.
#include "AssemblyGraph.hpp"
#include "Anchor.hpp"
#include "AnchorGraph.hpp"
#include "deduplicate.hpp"
#include "findLinearChains.hpp"
using namespace shasta;

// Standard library.
#include <fstream.hpp>
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

    cout << "The initial assembly graph has " << num_vertices(*this) <<
        " vertices and " << num_edges(*this) << " edges." << endl;
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
