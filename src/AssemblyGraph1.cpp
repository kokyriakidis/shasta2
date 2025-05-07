// Sjasta2.
#include "AssemblyGraph1.hpp"
#include "findLinearChains.hpp"
#include "TransitionGraph.hpp"
using namespace shasta;

// Standard library.
#include <fstream.hpp>
#include <map>


AssemblyGraph1::AssemblyGraph1(
    const AnchorGraph&,
    const TransitionGraph&)
{
    SHASTA_ASSERT(0);   // Use AssemblyGraph2 instead.

#if 0
    AssemblyGraph1& assemblyGraph1 = *this;


    // Find linear chains in the TransitionGraph.
    vector< vector<TransitionGraph::edge_descriptor> > chains;
    findLinearChains(transitionGraph, 0, chains);
    cout << "Found " << chains.size() << " linear chains in the TransitionGraph." << endl;


    // The initial and final vertices of each linear chain generate AssemblyGraph1 vertices.
    std::map<TransitionGraph::vertex_descriptor, vertex_descriptor> vertexMap;
    for(const vector<TransitionGraph::edge_descriptor>& chain: chains) {

        const TransitionGraph::edge_descriptor e0 = chain.front();
        const TransitionGraph::edge_descriptor e1 = chain.back();

        const TransitionGraph::vertex_descriptor v0 = source(e0, transitionGraph);
        const TransitionGraph::vertex_descriptor v1 = target(e1, transitionGraph);

        // Get the source vertex, creating it if necessary.
        vertex_descriptor u0 = null_vertex();
        auto it0 = vertexMap.find(v0);
        if(it0 == vertexMap.end()) {
            u0 = add_vertex(AssemblyGraph1Vertex(nextVertexId++), assemblyGraph1);
            vertexMap.insert(make_pair(v0, u0));
        } else {
            u0 = it0->second;
        }

        // Get the target vertex, creating it if necessary.
        vertex_descriptor u1 = null_vertex();
        auto it1 = vertexMap.find(v1);
        if(it1 == vertexMap.end()) {
            u1 = add_vertex(AssemblyGraph1Vertex(nextVertexId++),  assemblyGraph1);
            vertexMap.insert(make_pair(v1, u1));
        } else {
            u1 = it1->second;
        }

        // Create the edge.
        edge_descriptor e;
        bool edgeWasAdded;
        tie(e, edgeWasAdded) = add_edge(u0, u1, AssemblyGraph1Edge(nextEdgeId++), assemblyGraph1);
        SHASTA_ASSERT(edgeWasAdded);
        AssemblyGraph1Edge& edge = assemblyGraph1[e];
        for(const TransitionGraph::edge_descriptor eTransition: chain) {
            const TransitionGraph::vertex_descriptor vTransition = source(eTransition, transitionGraph);
            const AnchorGraph::edge_descriptor eAnchorGraph = transitionGraph[vTransition].eAnchorGraph;
            const AnchorPair& anchorPair = anchorGraph[eAnchorGraph];
            edge.push_back(anchorPair);
        }
        const AnchorGraph::edge_descriptor eAnchorGraph = transitionGraph[v1].eAnchorGraph;
        const AnchorPair& anchorPair = anchorGraph[eAnchorGraph];
        edge.push_back(anchorPair);
    }

    cout << "The initial AssemblyGraph1 has " << num_vertices(assemblyGraph1) <<
        " vertices and " << num_edges(assemblyGraph1) << " edges." << endl;
#endif
}



void AssemblyGraph1::writeGfa(const string& fileName) const
{
    ofstream gfa(fileName);
    writeGfa(gfa);
}



void AssemblyGraph1::writeGfa(ostream& gfa) const
{
    const AssemblyGraph1& assemblyGraph1 = *this;

    // Write the header line.
    gfa << "H\tVN:Z:1.0\n";

    // Write a segment for each edge.
    BGL_FORALL_EDGES(e, assemblyGraph1, AssemblyGraph1) {
        const AssemblyGraph1Edge& edge = assemblyGraph1[e];

        // Record type.
        gfa << "S\t";

        // Name.
        gfa << edge.edgeId << "\t";

        // Sequence.
        gfa << "*\t";
        gfa << "LN:i:" << 1000 * edge.size() << "\n";   // For now
    }



    // For each vertex, write links between each pair of incoming/outgoing edges.
    BGL_FORALL_VERTICES(v, assemblyGraph1, AssemblyGraph1) {
        BGL_FORALL_INEDGES(v, e0, assemblyGraph1, AssemblyGraph1) {
            const uint64_t id0 = assemblyGraph1[e0].edgeId;
            BGL_FORALL_OUTEDGES(v, e1, assemblyGraph1, AssemblyGraph1) {
                const uint64_t id1 = assemblyGraph1[e1].edgeId;
                gfa <<
                    "L\t" <<
                    id0 << "\t+\t" <<
                    id1 << "\t+\t*\n";
            }
        }
    }

}
