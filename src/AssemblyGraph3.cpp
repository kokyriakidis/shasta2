// Shasta.
#include "AssemblyGraph3.hpp"
#include "Anchor.hpp"
#include "AnchorGraph.hpp"
#include "findLinearChains.hpp"
using namespace shasta;

// Standard library.
#include <fstream.hpp>

// Explicit instantiation.
#include "MultithreadedObject.tpp"
template class MultithreadedObject<AssemblyGraph3>;



// Initial construction from the AnchorGraph.
AssemblyGraph3::AssemblyGraph3(
    const Anchors& anchors,
    const AnchorGraph& anchorGraph,
    const AssemblerOptions& assemblerOptions) :
    MappedMemoryOwner(anchors),
    MultithreadedObject<AssemblyGraph3>(*this),
    anchors(anchors),
    assemblerOptions(assemblerOptions)
{
    AssemblyGraph3& assemblyGraph3 = *this;

    // Find linear chains of edges in the AnchorGraph.
    vector< std::list<AnchorGraph::edge_descriptor> > chains;
    findLinearChains(anchorGraph, 1, chains);

    // Generate vertices.
    // At this stage there is a vertex for each AnchorGraph vertex
    // that is at the beginning or end of a linear chain,
    // so there is only one vertex for a given AnchorId.
    // However, after detangling there can be more than one vertex
    // with a given AnchorId. So the vertexMap is only used in this constructor.
    std::map<AnchorId, vertex_descriptor> vertexMap;
    for(const auto& chain: chains) {
        const AnchorId anchorId0 = anchorGraph[chain.front()].anchorPair.anchorIdA;
        const AnchorId anchorId1 = anchorGraph[chain.back()].anchorPair.anchorIdB;

        if(not vertexMap.contains(anchorId0)) {
            const vertex_descriptor v0 = add_vertex(AssemblyGraph3Vertex(anchorId0, nextVertexId++), assemblyGraph3);
            vertexMap.insert(make_pair(anchorId0, v0));
        }

        if(not vertexMap.contains(anchorId1)) {
            const vertex_descriptor v1 = add_vertex(AssemblyGraph3Vertex(anchorId1, nextVertexId++), assemblyGraph3);
            vertexMap.insert(make_pair(anchorId1, v1));
        }
    }
    SHASTA_ASSERT(vertexMap.size() == num_vertices(assemblyGraph3));



    // Generate the edges. There is an edge for each linear chain.
    for(const auto& chain: chains) {
        const AnchorId anchorId0 = anchorGraph[chain.front()].anchorPair.anchorIdA;
        const AnchorId anchorId1 = anchorGraph[chain.back()].anchorPair.anchorIdB;

        const vertex_descriptor v0 = vertexMap.at(anchorId0);
        const vertex_descriptor v1 = vertexMap.at(anchorId1);

        edge_descriptor e;
        bool edgeWasAdded;
        tie(e, edgeWasAdded) = add_edge(v0, v1, AssemblyGraph3Edge(nextEdgeId++), assemblyGraph3);
        AssemblyGraph3Edge& edge = assemblyGraph3[e];

        // Each AnchorGraph edge in the chain contributes a step to this AssemblyGraph3 edge.
        for(const AnchorGraph::edge_descriptor eA: chain) {
            const AnchorGraphEdge& edgeA = anchorGraph[eA];
            edge.emplace_back(edgeA.anchorPair, edgeA.offset);
        }
    }


    cout << "The AssemblyGraph3 has " << num_vertices(assemblyGraph3) <<
        " vertices and " << num_edges(assemblyGraph3) << " edges." << endl;

    writeGfa("AssemblyGraph3.gfa");
    check();
}



void AssemblyGraph3::check() const
{
    const AssemblyGraph3& assemblyGraph3 = *this;

    BGL_FORALL_EDGES(e, assemblyGraph3, AssemblyGraph3) {
        const AssemblyGraph3Edge& edge = assemblyGraph3[e];
        SHASTA_ASSERT(not edge.empty());



        // Check that the first/last AnchorIds of this edge are consistent
        // with the ones in the source/target vertices.
        const vertex_descriptor v0 = source(e, assemblyGraph3);
        const vertex_descriptor v1 = target(e, assemblyGraph3);

        const AnchorId anchorId0 = assemblyGraph3[v0].anchorId;
        const AnchorId anchorId1 = assemblyGraph3[v1].anchorId;

        SHASTA_ASSERT(edge.front().anchorPair.anchorIdA == anchorId0);
        SHASTA_ASSERT(edge.back().anchorPair.anchorIdB == anchorId1);



        // Check that AnchorPairs in this edge are adjacent to each other.
        for(uint64_t i1=1; i1<edge.size(); i1++) {
            const uint64_t i0 = i1 - 1;
            SHASTA_ASSERT(edge[i0].anchorPair.anchorIdB == edge[i1].anchorPair.anchorIdA);
        }    }

}



uint64_t AssemblyGraph3Edge::offset() const
{
    uint64_t sum = 0;
    for(const auto& step: *this) {
        sum += step.offset;
    }
    return sum;
}



void AssemblyGraph3::writeGfa(const string& fileName) const
{
    const AssemblyGraph3& assemblyGraph3 = *this;
    ofstream gfa(fileName);

    // Write the header line.
    gfa << "H\tVN:Z:1.0\n";

    // Each edge generates a gfa segment.
    BGL_FORALL_EDGES(e, assemblyGraph3, AssemblyGraph3) {
        const AssemblyGraph3Edge& edge = assemblyGraph3[e];

        // Record type.
        gfa << "S\t";

        // Name.
        gfa << edge.id << "\t";

        // Sequence.
        gfa << "*\t";
        gfa << "LN:i:" << edge.offset() << "\n";
    }



    // For each vertex, generate a link between each pair of
    // incoming/outgoing edges.
    BGL_FORALL_VERTICES(v, assemblyGraph3, AssemblyGraph3) {
        BGL_FORALL_INEDGES(v, e0, assemblyGraph3, AssemblyGraph3) {
            const uint64_t id0 = assemblyGraph3[e0].id;
            BGL_FORALL_OUTEDGES(v, e1, assemblyGraph3, AssemblyGraph3) {
                const uint64_t id1 = assemblyGraph3[e1].id;

                gfa <<
                    "L\t" <<
                    id0 << "\t+\t" <<
                    id1 << "\t+\t*\n";
            }
        }
    }


}
