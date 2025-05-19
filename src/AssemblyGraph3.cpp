// Shasta.
#include "AssemblyGraph3.hpp"
#include "Anchor.hpp"
#include "AnchorGraph.hpp"
#include "AssemblerOptions.hpp"
#include "findLinearChains.hpp"
#include "LocalAssembly2.hpp"
#include "performanceLog.hpp"
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


    check();
}



// Detangle, phase, assemble sequence, output.
void AssemblyGraph3::run(uint64_t threadCount)
{
    AssemblyGraph3& assemblyGraph3 = *this;

    cout << "The AssemblyGraph3 has " << num_vertices(assemblyGraph3) <<
        " vertices and " << num_edges(assemblyGraph3) << " edges." << endl;

    writeGfa("AssemblyGraph3-A.gfa");

    assembleAll(threadCount);
    writeGfa("AssemblyGraph3-Z.gfa");
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
    vector<shasta::Base> sequence;
    BGL_FORALL_EDGES(e, assemblyGraph3, AssemblyGraph3) {
        const AssemblyGraph3Edge& edge = assemblyGraph3[e];

        // Record type.
        gfa << "S\t";

        // Name.
        gfa << edge.id << "\t";

        // Sequence.
        if(edge.wasAssembled) {
            edge.getSequence(sequence);
            copy(sequence.begin(), sequence.end(), ostream_iterator<shasta::Base>(gfa));
            gfa << "\t";
            gfa << "LN:i:" << sequence.size() << "\n";

        } else {
            gfa << "*\t";
            gfa << "LN:i:" << edge.offset() << "\n";
        }
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



// Assemble sequence for all edges.
void AssemblyGraph3::assembleAll(uint64_t threadCount)
{
    performanceLog << timestamp << "Sequence assembly begins." << endl;

    const AssemblyGraph3& assemblyGraph3 = *this;

    edgesToBeAssembled.clear();
    BGL_FORALL_EDGES(e, assemblyGraph3, AssemblyGraph3) {
        edgesToBeAssembled.push_back(e);
    }
    assemble(threadCount);
    edgesToBeAssembled.clear();

    performanceLog << timestamp << "Sequence assembly ends." << endl;
}



// Assemble sequence for the specified edge.
void AssemblyGraph3::assemble(edge_descriptor e, uint64_t threadCount)
{
    edgesToBeAssembled.clear();
    edgesToBeAssembled.push_back(e);
    assemble(threadCount);
}



// Assemble sequence for step i of the specified edge.
// This is the lowest level sequence assembly function and is not multithreaded.
// It runs a LocalAssembly2 on the AnchorPair for that step.
void AssemblyGraph3::assembleStep(edge_descriptor e, uint64_t i)
{
    AssemblyGraph3& assemblyGraph3 = *this;
    AssemblyGraph3Edge& edge = assemblyGraph3[e];
    AssemblyGraph3EdgeStep& step = edge[i];

    // Run the LocalAssembly2.
    ofstream html;  // Not open, so no html output takes place.
    LocalAssembly2 localAssembly(
        anchors, html, false,
        assemblerOptions.aDrift,
        assemblerOptions.bDrift,
        step.anchorPair);
    localAssembly.run(false, assemblerOptions.localAssemblyOptions.maxAbpoaLength);
    localAssembly.getSequence(step.sequence);
}



// Assemble sequence for all edges in the edgesToBeAssembled vector.
// This fills in the stepsToBeAssembled with all steps of those edges,
// then assembles each of the steps in parallel.
void AssemblyGraph3::assemble(uint64_t threadCount)
{
    AssemblyGraph3& assemblyGraph3 = *this;

    stepsToBeAssembled.clear();
    for(const edge_descriptor e: edgesToBeAssembled) {
        AssemblyGraph3Edge& edge = assemblyGraph3[e];
        for(uint64_t i=0; i<edge.size(); i++) {
            stepsToBeAssembled.push_back(make_pair(e, i));
        }
    }

    const uint64_t batchCount = 1;
    setupLoadBalancing(stepsToBeAssembled.size(), batchCount);
    runThreads(&AssemblyGraph3::assembleThreadFunction, threadCount);

    // Mark them as assembled.
    for(const edge_descriptor e: edgesToBeAssembled) {
        assemblyGraph3[e].wasAssembled = true;
    }

    edgesToBeAssembled.clear();
    stepsToBeAssembled.clear();
}



void AssemblyGraph3::assembleThreadFunction(uint64_t /* threadId */)
{
    AssemblyGraph3& assemblyGraph3 = *this;

    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over all assembly steps assigned to this batch.
        for(uint64_t j=begin; j!=end; j++) {
            if((j % 1000) == 0) {
                std::lock_guard<std::mutex> lock(mutex);
                performanceLog << timestamp << j << "/" << stepsToBeAssembled.size() << endl;
            }

            const auto& p = stepsToBeAssembled[j];
            const edge_descriptor e = p.first;
            const uint64_t i = p.second;
            AssemblyGraph3Edge& edge = assemblyGraph3[e];
            SHASTA_ASSERT(i < edge.size());
            assembleStep(e, i);
        }
    }
}



// Clear sequence from all steps of all edges.
void AssemblyGraph3::clearSequence()
{
    AssemblyGraph3& assemblyGraph3 = *this;

    BGL_FORALL_EDGES(e, assemblyGraph3, AssemblyGraph3) {
        AssemblyGraph3Edge& edge = assemblyGraph3[e];
        edge.wasAssembled= false;
        for(AssemblyGraph3EdgeStep& step: edge) {
            step.sequence.clear();
            step.sequence.shrink_to_fit();
        }
    }
}



void AssemblyGraph3Edge::getSequence(vector<Base>& sequence) const
{
    sequence.clear();
    for(const auto& step: *this) {
        copy(step.sequence.begin(), step.sequence.end(), back_inserter(sequence));
    }
}
