#pragma once

// In the AssemblyGraph, each edge is a linear chains of adjacent AnchorPairs.
// It is initially created from linear chains of edges in the AnchorGraph.
// Each AssemblyGraphEdge generates a gfa segment.

// Shasta.
#include "AnchorPair.hpp"
#include "Base.hpp"
#include "invalid.hpp"
#include "MappedMemoryOwner.hpp"
#include "MultithreadedObject.hpp"

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>
#include <boost/serialization/vector.hpp>

// Standard library.
#include <tuple.hpp>


namespace shasta {

    class AssemblyGraph;
    class AssemblyGraphVertex;
    class AssemblyGraphEdge;
    class AssemblyGraphEdgeStep;

    using AssemblyGraphBaseClass = boost::adjacency_list<
        boost::listS,
        boost::listS,
        boost::bidirectionalS,
        AssemblyGraphVertex,
        AssemblyGraphEdge>;

    class AnchorGraph;
    class Anchors;
    class Options;
    class Detangler;
    class Superbubble;
    class SuperbubbleChain;

}



// Tha AnchorId of a vertex is the last AnchorId of the last AnchorPair
// of all incoming edges and the first AnchorId of the first AnchorPair
// of each outgoing edge.
// When the AssemblyGraph is initially created from the AnchorGraph,
// there can be at most one vertex for each AnchorId.
// However that is no longer true after detangling.
class shasta::AssemblyGraphVertex {
public:
    AnchorId anchorId = invalid<AnchorId>;
    uint64_t id = invalid<uint64_t>;

    AssemblyGraphVertex(AnchorId anchorId, uint64_t id) :
        anchorId(anchorId), id(id) {}
    AssemblyGraphVertex() {}

    template<class Archive> void serialize(Archive& ar, unsigned int /* version */)
    {
        ar & anchorId;
        ar & id;
    }
};



class shasta::AssemblyGraphEdgeStep {
public:
    AnchorPair anchorPair;
    uint64_t offset = invalid<uint64_t>;
    vector<Base> sequence;

    AssemblyGraphEdgeStep(
        const AnchorPair& anchorPair,
        uint64_t offset) :
        anchorPair(anchorPair),
        offset(offset)
    {}

    AssemblyGraphEdgeStep()
    {}

    template<class Archive> void serialize(Archive& ar, unsigned int /* version */)
    {
        ar & anchorPair;
        ar & offset;
        ar & sequence;
    }
};



class shasta::AssemblyGraphEdge : public vector<AssemblyGraphEdgeStep> {
public:
    uint64_t id = invalid<uint64_t>;
    bool wasAssembled = false;

    AssemblyGraphEdge(uint64_t id = invalid<uint64_t>) : id(id) {}

    void check(const Anchors&) const;

    uint64_t offset() const;
    uint64_t sequenceLength() const;
    void getSequence(vector<Base>&) const;

    double averageCoverage() const;

    template<class Archive> void serialize(Archive& ar, unsigned int /* version */)
    {
        ar & boost::serialization::base_object< vector<AssemblyGraphEdgeStep> >(*this);
        ar & id;
        ar & wasAssembled;
    }

    void swapSteps(AssemblyGraphEdge& that)
    {
        vector<AssemblyGraphEdgeStep>::swap(that);
    }
};



class shasta::AssemblyGraph :
    public AssemblyGraphBaseClass,
    public MappedMemoryOwner,
    public MultithreadedObject<AssemblyGraph> {
public:

    // Initial construction from the AnchorGraph.
    AssemblyGraph(
        const Anchors&,
        const Journeys&,
        const AnchorGraph&,
        const Options&);

    // Deserialize constructor.
    AssemblyGraph(
        const Anchors&,
        const Journeys&,
        const Options&,
        const string& stage);

    // Detangle, phase, assemble sequence, output.
    void run(uint64_t threadCount);

    const Anchors& anchors;
    const Journeys& journeys;
    uint64_t nextVertexId = 0;
    uint64_t nextEdgeId = 0;
    const Options& options;

private:

    void check() const;

    // Bubble cleanup.
    // A bubble is a set of parallel edges in the AssemblyGraph.
    class Bubble {
    public:
        vertex_descriptor v0;
        vertex_descriptor v1;
        vector<edge_descriptor> edges;
    };
    void findBubbles(vector<Bubble>&) const;
    void bubbleCleanup(uint64_t threadCount);
    uint64_t bubbleCleanupIteration(uint64_t threadCount);

    void prune();

public:
    // Compress linear chains of edges into a single edge.
    void compress();

    // Class to order vertices or edges by id.
    class OrderById {
    public:
        OrderById(const AssemblyGraph& assemblyGraph): assemblyGraph(assemblyGraph) {}
        const AssemblyGraph& assemblyGraph;
        bool operator()(vertex_descriptor x, vertex_descriptor y) const
        {
            return assemblyGraph[x].id < assemblyGraph[y].id;
        }
        bool operator()(edge_descriptor x, edge_descriptor y) const
        {
            return assemblyGraph[x].id < assemblyGraph[y].id;
        }
    };

    bool hasSelfEdge(vertex_descriptor v) const
    {
        bool edgeExists = false;
        tie(ignore, edgeExists) = boost::edge(v, v, *this);
        return edgeExists;
    }


    // Detangling.
    // The detangling functions return the number of successful detangling operations.

    // Low level detangling function.
    uint64_t detangle(const vector< vector<vertex_descriptor> >& detanglingCandidates, Detangler&);

    uint64_t detangleVerticesIteration(Detangler&);
    uint64_t detangleVertices(uint64_t maxIterationCount, Detangler&);

    uint64_t detangleEdgesIteration(
        uint64_t maxEdgeLength,
        Detangler&);
    uint64_t detangleEdges(
        uint64_t maxIterationCount,
        uint64_t maxEdgeLength,
        Detangler&);

    uint64_t detangleTemplateIteration(uint64_t templateId, Detangler&);
    uint64_t detangleTemplate(uint64_t templateId, uint64_t maxIterationCount, Detangler&);
    uint64_t detangleTemplates(uint64_t maxIterationCount, Detangler&);

    // High level detangling function.
    uint64_t detangle(
        uint64_t maxIterationCount,
        uint64_t maxEdgeLength,
        Detangler&);



    // Tangle templates.
    using TangleTemplate = boost::adjacency_list<boost::listS, boost::vecS, boost::bidirectionalS>;
    static void writeGraphviz(ostream&, const TangleTemplate&);
    static TangleTemplate reverse(const TangleTemplate&);
    vector<TangleTemplate> tangleTemplates;
    void createTangleTemplates();


public:

    // Compute oriented read journeys in the AssemblyGraph.
    void computeJourneys() const;


    // Superbubbles and SuperbubbleChains.
    // - A Superbubble is a Tangle in which all entrance edges are into
    //   a single vertex and all exit edges are from a single vertex.
    // - A SuperbubbleChain is a linear sequence of Superbubbles.
    //   In the sequence, the target vertex of a Superbubble is the same
    //   as the source vertex of the next Superbubble in the SuperbubbleChain.

    // This finds all Superbubbles seen using options.findSuperbubblesMaxDistance.
    // Some pairs of Superbubble can intersect (that is, they can have common edges).
    void findSuperbubbles(vector<Superbubble>&) const;

    // Remove Superbubbles that are entirely contained in a larger superbubble.
    void removeContainedSuperbubbles(vector<Superbubble>&) const;

    void writeSuperbubbles(const vector<Superbubble>&, const string& fileName) const;
    void writeSuperbubblesForBandage(const vector<Superbubble>&, const string& fileName) const;

    void findSuperbubbleChains(
        const vector<Superbubble>&,
        vector<SuperbubbleChain>&
        ) const;

    void writeSuperbubbleChains(const vector<SuperbubbleChain>&, const string& fileName) const;
    void writeSuperbubbleChainsForBandage(const vector<SuperbubbleChain>&, const string& fileName) const;

    // Simplify Superbubbles by turning them into bubbles via clustering
    // of oriented read journeys.
    void simplifySuperbubbles();
    void simplifySuperbubble(const Superbubble&, uint64_t minCoverage);

    // Phasing of SuperbubbleChains.
    void phaseSuperbubbleChains();



    // Strongly connected components.

    // Find the non-trivial strongly connected components.
    // Each component is stored with vertices sorted to permit binary searches.
    void findStrongComponents(vector< vector<vertex_descriptor> >&) const;

    // This creates a csv file that can be loaded in bandage to see
    // the strongly connected components.
    void colorStrongComponents() const;


    // Output.
    void write(const string& stage);
private:
    void writeGfa(const string& fileName) const;
    void writeGfa(ostream&) const;
    void writeFasta(const string& stage) const;
    void writeGraphviz(const string& fileName) const;
    void writeGraphviz(ostream&) const;



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
        ar & boost::serialization::base_object<AssemblyGraphBaseClass>(*this);
        ar & nextVertexId;
        ar & nextEdgeId;
    }
    void save(ostream&) const;
    void load(istream&);

    // These do save/load to/from mapped memory.
    // The file name is AssemblyGraph-Stage.
    void save(const string& stage) const;
    void load(const string& stage);
};
