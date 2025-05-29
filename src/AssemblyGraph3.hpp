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
#include <boost/serialization/vector.hpp>


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
    class Detangler;

}



// Tha AnchorId of a vertex is the last AnchorId of the last AnchorPair
// of all incoming edges and the first AnchorId of the first AnchorPair
// of each outgoing edge.
// When the Assembly3Graph is initially created from the AnchorGraph,
// there can be at most one vertex for each AnchorId.
// However that is no longer true after detangling.
class shasta::AssemblyGraph3Vertex {
public:
    AnchorId anchorId = invalid<AnchorId>;
    uint64_t id = invalid<uint64_t>;

    AssemblyGraph3Vertex(AnchorId anchorId, uint64_t id) :
        anchorId(anchorId), id(id) {}
    AssemblyGraph3Vertex() {}

    template<class Archive> void serialize(Archive& ar, unsigned int /* version */)
    {
        ar & anchorId;
        ar & id;
    }
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

    template<class Archive> void serialize(Archive& ar, unsigned int /* version */)
    {
        ar & anchorPair;
        ar & offset;
        ar & sequence;
    }
};



class shasta::AssemblyGraph3Edge : public vector<AssemblyGraph3EdgeStep> {
public:
    uint64_t id = invalid<uint64_t>;
    bool wasAssembled = false;

    AssemblyGraph3Edge(uint64_t id = invalid<uint64_t>) : id(id) {}

    void check(const Anchors&) const;

    uint64_t offset() const;
    uint64_t sequenceLength() const;
    void getSequence(vector<Base>&) const;

    template<class Archive> void serialize(Archive& ar, unsigned int /* version */)
    {
        ar & boost::serialization::base_object< vector<AssemblyGraph3EdgeStep> >(*this);
        ar & id;
        ar & wasAssembled;
    }

    void swapSteps(AssemblyGraph3Edge& that)
    {
        vector<AssemblyGraph3EdgeStep>::swap(that);
    }
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

    // Deserialize constructor.
    AssemblyGraph3(
        const Anchors&,
        const AssemblerOptions&,
        const string& stage);

    // Detangle, phase, assemble sequence, output.
    void run(uint64_t threadCount);

    const Anchors& anchors;
    uint64_t nextVertexId = 0;
    uint64_t nextEdgeId = 0;
private:
    const AssemblerOptions& assemblerOptions;

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

public:
    // Compress linear chains of edges into a single edge.
    void compress();



    // Detangling.

    // One iteration of detangling a given set of tangles using the given Detangler.
    // This is the lowest level detangling function.
    bool detangle(const vector< vector<vertex_descriptor> >& detanglingCandidates, Detangler&);

    // One iteration of vertex detangling using the given Detangler.
    bool detangleVertices(Detangler&);

    // One iteration of edge detangling using the given Detangler.
    bool detangleEdges(Detangler&);

    // One iteration of all usable detangling functions using the given Detangler.
    bool detangleIteration(Detangler&);

    // Multiple iterations of all usable detangling functions using the given Detangler.
    bool detangle(uint64_t maxIterationCount, Detangler&);



    // Output.
public:
    void write(const string& stage);
private:
    void writeGfa(const string& fileName) const;
    void writeGfa(ostream&) const;
    void writeFasta(const string& stage) const;



    // Sequence assembly.

    // Assemble sequence for all edges.
public:
    void assembleAll(uint64_t threadCount);
private:

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



    // Serialization.
    friend class boost::serialization::access;
    template<class Archive> void serialize(Archive& ar, unsigned int /* version */)
    {
        ar & boost::serialization::base_object<AssemblyGraph3BaseClass>(*this);
        ar & nextVertexId;
        ar & nextEdgeId;
    }
    void save(ostream&) const;
    void load(istream&);

    // These do save/load to/from mapped memory.
    // The file name is AssemblyGraph3-Stage.
    void save(const string& stage) const;
    void load(const string& stage);
};
