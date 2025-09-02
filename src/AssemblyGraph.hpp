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
#include <boost/graph/adj_list_serialize.hpp>
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



// The AnchorId of a vertex is the last AnchorId of the last AnchorPair
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

    // A call to countOrientedReadStepsBySegment will store here the distinct OrientedReadIds
    // that visit this edge and that also visit at least one other edge.
    // They are stored sorted, each with the number of steps in this edge
    // in which the OrientedReadId appears.
    // This is used to compute extended tangle matrices.
    vector< pair<OrientedReadId, uint64_t> > transitioningOrientedReadIds;

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

    AnchorId firstAnchorId() const
    {
        return front().anchorPair.anchorIdA;
    }
    AnchorId lastAnchorId() const
    {
        return back().anchorPair.anchorIdB;
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
    void simplifyAndAssemble();

    const Anchors& anchors;
    const Journeys& journeys;
    uint64_t nextVertexId = 0;
    uint64_t nextEdgeId = 0;
    const Options& options;

private:

    void check() const;

    // Count the distinct target vertices among
    // all the out-edges of a given vertex.
    // So this counts the number of parallel edge sets
    // outgoing from the given vertex.
    uint64_t countDistinctTargetVertices(vertex_descriptor) const;

    // Count the distinct source vertices among
    // all the in-edges of a given vertex.
    // So this counts the number of parallel edge sets
    // incoming into the given vertex.
    uint64_t countDistinctSourceVertices(vertex_descriptor) const;

    // Bubble cleanup.
    // A bubble is a set of parallel edges in the AssemblyGraph.
    // The edges are sorted by id.
    class Bubble {
    public:
        vertex_descriptor v0;
        vertex_descriptor v1;
        vector<edge_descriptor> edges;
    };
    void findBubbles(vector<Bubble>&) const;
public:
    uint64_t bubbleCleanup();
    uint64_t bubbleCleanupIteration(vector< pair<vertex_descriptor, vertex_descriptor> >& excludeList);
private:
    bool bubbleCleanup(const Bubble&);


    // Analyze a Bubble and finds pairs of "similar" branches.
    // See analyzeSimilarSequences.h for more information.
    bool analyzeBubble(
        const Bubble&,
        const vector<uint64_t> minRepeatCount,
        vector< pair<uint64_t, uint64_t> >& similarPairs
        ) const;

public:

    uint64_t prune();

    // Compress linear chains of edges into a single edge.
    uint64_t compress();
    uint64_t compressDebugLevel = 0; // 1=compact, 2=detailed.



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

        // Also order pairs of edges.
        using EdgePair = pair<edge_descriptor, edge_descriptor>;
        bool operator()(const EdgePair& x,const EdgePair& y) const
        {
            if(assemblyGraph[x.first].id < assemblyGraph[y.first].id) {
                return true;
            }
            if(assemblyGraph[x.first].id > assemblyGraph[y.first].id) {
                return false;
            }
            return assemblyGraph[x.second].id < assemblyGraph[y.second].id;
        }
    };
    const OrderById orderById;



    bool hasSelfEdge(vertex_descriptor v) const
    {
        bool edgeExists = false;
        tie(ignore, edgeExists) = boost::edge(v, v, *this);
        return edgeExists;
    }



    // Detangling.
    // The detangling functions return the number of successful detangling operations.

    // High level detangling function.
    // The second detangler is used for short tangles.
    uint64_t detangleHighLevel(Detangler&, Detangler&);

    // Vertex detangling.
    uint64_t detangleVerticesIteration(Detangler&);
    uint64_t detangleVertices(Detangler&);

    // Edge detangling.
    uint64_t detangleEdgesIteration(Detangler&);
    uint64_t detangleEdges(Detangler&);

    // Detangling of short tangles.
    uint64_t detangleShortTangles(Detangler&);

    // Low level detangling function.
    uint64_t detangleLowLevel(const vector< vector<vertex_descriptor> >& detanglingCandidates, Detangler&);


    // The detangling process can generate empty edges (edges without steps).
    // This removes them by collapsing the vertices they join.
    void removeEmptyEdges();



public:

    // Compute compressed journeys in the AssemblyGraph.
    // The compressed journey of an oriented read
    // is the sequence of assembly graph edges it visits.
    void computeJourneys();
    vector< vector<edge_descriptor> > compressedJourneys;



    // Find appearances of each OrientedReadId in the AssemblyGraphSteps
    // of each AssemblyGraphEdge.
    // This stores information for each OrientedReadId in
    // orientedReadEdgeInformation (indexed by OrientedReadId::getValue).
    // This also stores in the each AssemblyGraphEdge::transitioningOrientedReadIds
    // the OrientedReadIds that visit the edge and at least one other edge.
    void findOrientedReadEdgeInformation();
    void clearOrientedReadEdgeInformation();
    void writeOrientedReadEdgeInformation();
    class OrientedReadEdgeInformation {
    public:
        edge_descriptor e;
        uint64_t stepCount;
        OrientedReadEdgeInformation(edge_descriptor e, uint64_t stepCount) :
            e(e), stepCount(stepCount) {}
    };
    class OrientedReadEdgeInformationOrderById {
    public:
        OrientedReadEdgeInformationOrderById(const AssemblyGraph& assemblyGraph): assemblyGraph(assemblyGraph) {}
        const AssemblyGraph& assemblyGraph;
        bool operator()(const OrientedReadEdgeInformation& x, const OrientedReadEdgeInformation& y) const
        {
            return assemblyGraph.orderById(x.e, y.e);
        }
    };
    // Indexed by OrientedReadId::getValue()
    // For each OrientedReadId, the OrientedReadEdgeInformation are ordered by edge id.
    vector< vector<OrientedReadEdgeInformation> > orientedReadEdgeInformation;



    // Use the orientedReadSegments and the transitioningOrientedReadIds
    // stored in the AssemblyGraphEdges to compute an extended tangle matrix.
    void computeExtendedTangleMatrix(
        vector<edge_descriptor>& entrances,
        vector<edge_descriptor>& exits,
        vector< vector<double> > & tangleMatrix,
        ostream& html
        ) const;



    // Gather information on the oriented reads that appear
    // in the AssemblyGraphSteps of an edge.
    class OrientedReadEdgeInfo {
    public:
        OrientedReadId orientedReadId;
        uint64_t stepCount;
        uint32_t minPositionInJourney;
        uint32_t maxPositionInJourney;
        bool operator<(const OrientedReadEdgeInfo& that) const
        {
            return orientedReadId < that.orientedReadId;
        }
    };
    void gatherOrientedReadInformationOnEdge(
        edge_descriptor,
        vector<OrientedReadEdgeInfo>&) const;


    // Simple search starting at a given edge (segment) and moving in the specified direction.
    void search(
        edge_descriptor,
        uint64_t direction) const;
    class SearchGraphVertex {
    public:
        AssemblyGraph::edge_descriptor e;
        SearchGraphVertex(AssemblyGraph::edge_descriptor e) : e(e) {}
    };
    class SearchGraphEdge {
    public:
        uint64_t tangleMatrix;
    };
    class SearchGraph : public boost::adjacency_list<
        boost::listS,
        boost::listS,
        boost::bidirectionalS,
        SearchGraphVertex,
        SearchGraphEdge> {
        public:
    };



    // Local search that continues as long as we have exactly one way to move.
    void forwardLocalSearch(
        edge_descriptor,
        uint64_t lowCoverageThreshold,
        uint64_t highCoverageThreshold,
        vector<edge_descriptor>&
    ) const;
    void backwardLocalSearch(
        edge_descriptor,
        uint64_t lowCoverageThreshold,
        uint64_t highCoverageThreshold,
        vector<edge_descriptor>&
    ) const;
    void testLocalSearch(
        uint64_t id,
        uint64_t direction,
        uint64_t lowCoverageThreshold,
        uint64_t highCoverageThreshold
        ) const;
    void createSearchGraph(
        uint64_t lowCoverageThreshold,
        uint64_t highCoverageThreshold);



    // More systematic search functionality that uses an index.
    void findEdgePairs(uint64_t minCoverage);
    std::map<edge_descriptor, vector<edge_descriptor> > edgePairsBySource;
    std::map<edge_descriptor, vector<edge_descriptor> > edgePairsByTarget;
    void testSearch(uint64_t edgeId, uint64_t direction, uint64_t minCoverage) const;



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

    // Simplify Superbubbles.
    uint64_t simplifySuperbubbles(Detangler&);
    bool simplifySuperbubbleByDetangling(const Superbubble&, Detangler&);
    bool simplifySuperbubbleByClustering(const Superbubble&, uint64_t minCoverage, uint64_t maxOffset);

    // Phasing of SuperbubbleChains.
    uint64_t phaseSuperbubbleChains();



    // Strongly connected components.

    // Find the non-trivial strongly connected components.
    // Each component is stored with vertices sorted to permit binary searches.
    void findStrongComponents(vector< vector<vertex_descriptor> >&) const;

    // This creates a csv file that can be loaded in bandage to see
    // the strongly connected components.
    void colorStrongComponents() const;


    // Output.
    void write(const string& stage);
    void writeFasta(const string& stage) const;
private:
    void writeGfa(const string& fileName) const;
    void writeGfa(ostream&) const;
    void writeGraphviz(const string& fileName) const;
    void writeGraphviz(ostream&) const;
    void writeCsv(const string& fileName) const;
    void writeCsv(ostream&) const;



    // Sequence assembly.

    // Assemble sequence for all edges.
public:
    void assembleAll();
private:

    // Assemble sequence for the specified edge.
    void assemble(edge_descriptor);

    // Assemble sequence for step i of the specified edge.
    // This is the lowest level sequence assembly function and is not multithreaded.
    // It runs a LocalAssembly2 on the AnchorPair for that step.
    void assembleStep(edge_descriptor, uint64_t i);

    // Assemble sequence for all edges in the edgesToBeAssembled vector.
    // This fills in the stepsToBeAssembled with all steps of those edges,
    // then assembles each of the steps in parallel.
    void assemble();
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
