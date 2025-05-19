#pragma once

// In the AssemblyGraph3, each edge is a linear chains of adjacent AnchorPairs.
// It is initially created from linear chains of edges in the AnchorGraph.
// Each Assembly3GraphEdge generates a gfa segment.

// Shasta.
#include "AnchorPair.hpp"
#include "invalid.hpp"
#include "MappedMemoryOwner.hpp"
#include "MultithreadedObject.hpp"

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>


namespace shasta {

    class AssemblyGraph3;
    class AssemblyGraph3Vertex;
    class AssemblyGraph3Edge;
    class AssemblyGraph3EdgeStep;

    using AssemblyGraph3BaseClass = boost::adjacency_list<
        boost::listS,
        boost::listS,
        boost::bidirectionalS,
        AssemblyGraph3Vertex,
        AssemblyGraph3Edge>;

    class AnchorGraph;
    class Anchors;
    class AssemblerOptions;

}



// Tha AnchorId of a vertex is the last AnchorId of the last AnchorPair
// of all incoming edges and the first AnchorId of the first AnchorPair
// of each outgoing edge.
class shasta::AssemblyGraph3Vertex {
public:
    AnchorId anchorId = invalid<AnchorId>;
    uint64_t id = invalid<uint64_t>;

    AssemblyGraph3Vertex(AnchorId anchorId, uint64_t id) :
        anchorId(anchorId), id(id) {}
};



class shasta::AssemblyGraph3EdgeStep {
public:
    AnchorPair anchorPair;
    uint64_t offset = invalid<uint64_t>;
    vector<Base> sequence;

    AssemblyGraph3EdgeStep(
        const AnchorPair& anchorPair,
        uint64_t offset) :
        anchorPair(anchorPair),
        offset(offset)
    {}

    AssemblyGraph3EdgeStep()
    {}

};



class shasta::AssemblyGraph3Edge : public vector<AssemblyGraph3EdgeStep> {
public:
    uint64_t id = invalid<uint64_t>;
    bool wasAssembled = false;

    AssemblyGraph3Edge(uint64_t id = invalid<uint64_t>) : id(id) {}

    void check(const Anchors&) const;

    uint64_t offset() const;
    void getSequence(vector<Base>&) const;
};



class shasta::AssemblyGraph3 :
    public AssemblyGraph3BaseClass,
    public MappedMemoryOwner,
    public MultithreadedObject<AssemblyGraph3> {
public:

    // Initial construction from the AnchorGraph.
    AssemblyGraph3(
        const Anchors&,
        const AnchorGraph&,
        const AssemblerOptions&);

    // Detangle, phase, assemble sequence, output.
    void run(uint64_t threadCount);

private:
    const Anchors& anchors;
    const AssemblerOptions& assemblerOptions;
    uint64_t nextVertexId = 0;
    uint64_t nextEdgeId = 0;

    void check() const;

    // Bubble cleanup.
    // A bubble is a set of parallel edges in the AssemblyGraph3.
    class Bubble {
    public:
        vertex_descriptor v0;
        vertex_descriptor v1;
        vector<edge_descriptor> edges;
    };
    void findBubbles(vector<Bubble>&) const;
    void bubbleCleanup(uint64_t threadCount);
    uint64_t bubbleCleanupIteration(uint64_t threadCount);

    // Compress linear chains of edges into a single edge.
    void compress();

    // Output.
    void write(const string& stage);
    void writeGfa(const string& fileName) const;
    void writeGfa(ostream&) const;
    void writeFasta(const string& stage) const;



    // Sequence assembly.

    // Assemble sequence for all edges.
    void assembleAll(uint64_t threadCount);

    // Assemble sequence for the specified edge.
    void assemble(edge_descriptor, uint64_t threadCount);

    // Assemble sequence for step i of the specified edge.
    // This is the lowest level sequence assembly function and is not multithreaded.
    // It runs a LocalAssembly2 on the AnchorPair for that step.
    void assembleStep(edge_descriptor, uint64_t i);

    // Assemble sequence for all edges in the edgesToBeAssembled vector.
    // This fills in the stepsToBeAssembled with all steps of those edges,
    // then assembles each of the steps in parallel.
    void assemble(uint64_t threadCount);
    void assembleThreadFunction(uint64_t threadId);
    vector<edge_descriptor> edgesToBeAssembled;
    vector< pair<edge_descriptor, uint64_t> > stepsToBeAssembled;

    // Clear sequence from all steps of all edges.
    void clearSequence();
};
