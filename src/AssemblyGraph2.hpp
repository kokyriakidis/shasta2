#pragma once

// In the AssemblyGraph2, each vertex is a linear chains of adjacent AnchorPairs.
// It is initially created from linear chains of vertices in the TransitionGraph.
// Each Assembly2GraphVertex generates a gfa segment.

// Shasta.
#include "AnchorPair.hpp"
#include "invalid.hpp"
#include "MappedMemoryOwner.hpp"
#include "MultithreadedObject.hpp"

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>
#include <boost/serialization/vector.hpp>

// Standard library.
#include <iosfwd.hpp>
#include <string.hpp>
#include <vector.hpp>

namespace shasta {

    class AssemblyGraph2;
    class AssemblyGraph2Vertex;
    class AssemblyGraph2VertexStep;
    class AssemblyGraph2Edge;

    using AssemblyGraph2BaseClass = boost::adjacency_list<
        boost::listS,
        boost::listS,
        boost::bidirectionalS,
        AssemblyGraph2Vertex,
        AssemblyGraph2Edge>;

    class Anchors;
    class AssemblerOptions;
    class Base;
    class TransitionGraph;
}



class shasta::AssemblyGraph2VertexStep {
public:
    AnchorPair anchorPair;
    uint64_t offset = invalid<uint64_t>;
    vector<Base> sequence;

    AssemblyGraph2VertexStep(
        const AnchorPair& anchorPair,
        uint64_t offset) :
        anchorPair(anchorPair),
        offset(offset)
    {}

    AssemblyGraph2VertexStep()
    {}

    template<class Archive> void serialize(Archive& ar, unsigned int /* version */)
    {
        ar & anchorPair;
        ar & offset;
        ar & sequence;
    }
};



class shasta::AssemblyGraph2Vertex : public vector<AssemblyGraph2VertexStep> {
public:
    uint64_t id = invalid<uint64_t>;
    bool wasAssembled = false;

    AssemblyGraph2Vertex(uint64_t id) : id(id) {}
    AssemblyGraph2Vertex() {}

    void check(const Anchors&) const;

    uint64_t offset() const;
    uint64_t sequenceLength() const;

    void getSequence(vector<Base>&) const;

    template<class Archive> void serialize(Archive& ar, unsigned int /* version */)
    {
        ar & boost::serialization::base_object< vector<AssemblyGraph2VertexStep> >(*this);
        ar & id;
        ar & wasAssembled;
    }
};



class shasta::AssemblyGraph2Edge {
public:
    template<class Archive> void serialize(Archive& /* ar */, unsigned int /* version */)
    {
    }
};



class shasta::AssemblyGraph2 :
    public AssemblyGraph2BaseClass,
    public MappedMemoryOwner,
    public MultithreadedObject<AssemblyGraph2> {
public:

    // Initial construction from the TransitionGraph.
    AssemblyGraph2(
        const Anchors&,
        const TransitionGraph&,
        const AssemblerOptions&);

    // Deserialize constructor.
    AssemblyGraph2(
        const Anchors&,
        const AssemblerOptions&,
        const string& stage);

    // Detangle, phase, assemble sequence, output.
    void run(uint64_t threadCount);

private:
    const Anchors& anchors;
    const AssemblerOptions& assemblerOptions;
    uint64_t nextVertexId = 0;

    void check() const;
    void check(edge_descriptor) const;

    // Bubble cleanup.
    // A bubble is a set of linear chains of vertices in the AssemblyGraph2.
    // In most cases each of the linear chains contains only one vertex.
    class Bubble {
    public:
        vertex_descriptor v0;
        vertex_descriptor v1;
        vector< vector<vertex_descriptor> > chains;
    };
    void findBubbles(vector<Bubble>&) const;
    void bubbleCleanup(uint64_t threadCount);

    // Output.
    void write(const string& stage);
    void writeGfa(const string& fileName) const;
    void writeGfa(ostream&) const;
    void writeFasta(const string& fileName) const;
    void writeGraphviz(const string& fileName) const;
    void writeGraphviz(ostream&) const;



    // Serialization.
    friend class boost::serialization::access;
    template<class Archive> void serialize(Archive& ar, unsigned int /* version */)
    {
        ar & boost::serialization::base_object<AssemblyGraph2BaseClass>(*this);
        ar & nextVertexId;
    }
    void save(ostream&) const;
    void load(istream&);

    // These do save/load to/from mapped memory.
    // The file name is AssemblyGraph2-Stage.
    void save(const string& stage) const;
    void load(const string& stage);



    // Merge vertices in linear chains.
    void compress();

    // Assemble sequence for all vertices.
    void assembleAll(uint64_t threadCount);

    // Assemble sequence for the specified vertex.
    void assemble(vertex_descriptor, uint64_t threadCount);

    // Assemble sequence for step i of the specified vertex.
    void assembleStep(vertex_descriptor, uint64_t i);

    // Assemble sequence for all vertices in the verticesToBeAssembled vector.
    // This fills in the stepToBeAssembled with all steps of those edges,
    // then assembles each of the steps in parallel.
    void assemble(uint64_t threadCount);
    void assembleThreadFunction(uint64_t threadId);
    vector<vertex_descriptor> verticesToBeAssembled;
    vector< pair<vertex_descriptor, uint64_t> > stepsToBeAssembled;

    void clearSequence();
};
