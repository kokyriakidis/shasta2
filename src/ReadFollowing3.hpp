#pragma once

// Read following in the AssemblyGraph.

// Shasta.
#include "AssemblyGraph.hpp"
#include "MultithreadedObject.hpp"
#include "SegmentStepSupport.hpp"

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>

// Standard library.
#include <map>
#include <random>
#include <set>



namespace shasta2 {

    namespace ReadFollowing3 {

        class Graph;
        class Vertex;
        class Edge;
        using GraphBaseClass = boost::adjacency_list<
            boost::listS,
            boost::listS,
            boost::bidirectionalS,
            Vertex,
            Edge>;

        using Path = vector<GraphBaseClass::vertex_descriptor>;

        // A Segment is an edge of the AssemblyGraph.
        using Segment = AssemblyGraph::edge_descriptor;
    }
}



class shasta2::ReadFollowing3::Vertex {
public:
    // A Segment is an edge of the AssemblyGraph.
    Segment segment;

    // The sequence length or estimated offset of this AssemblyGraph edge.
    uint64_t length = invalid<uint64_t>;

    // This is set for long vertices (length >= readFollowingSegmentLengthThreshold).
    bool isLong = false;

    Vertex(const AssemblyGraph&, Segment);
    Vertex() {}

    vector<SegmentStepSupport> initialSupport;
    vector<SegmentStepSupport> finalSupport;



    // Machinery to compute random paths biased by unnormalized edge probabilities.
    class MarkovTableItem {
    public:
        GraphBaseClass::vertex_descriptor v;
        double pUnnormalized = 0.;
        double p = 0.;
        double pCumulative = 0.;

        MarkovTableItem(
            GraphBaseClass::vertex_descriptor v,
            double pUnnormalized) :
            v(v), pUnnormalized(pUnnormalized) {}

        bool operator<(const MarkovTableItem& that) const
        {
            return p > that.p;
        }
    };
    array< vector<MarkovTableItem>, 2> markovTables;

    GraphBaseClass::vertex_descriptor next(uint64_t direction, double random01) const
    {
        for(const MarkovTableItem& markovTableItem: markovTables[direction]) {
            if(random01 <= markovTableItem.pCumulative) {
                return markovTableItem.v;
            }
        }
        SHASTA2_ASSERT(0);
    }



};



class shasta2::ReadFollowing3::Edge {
public:
    Edge(const AssemblyGraph&, Segment, Segment);
    SegmentPairInformation segmentPairInformation;

    // Unnormalized probability of this edge, for random paths.
    double pUnnormalized = 0.;
};



class shasta2::ReadFollowing3::Graph :
    public GraphBaseClass,
    public MultithreadedObject<Graph> {
public:
    Graph(const AssemblyGraph&);

    const AssemblyGraph& assemblyGraph;

    // Vertex creation.
    std::map<Segment, vertex_descriptor> vertexMap;
    void createVertices();

    // Create edge candidates.
    class EdgeCandidate {
    public:
        vertex_descriptor v0;
        vertex_descriptor v1;
        bool operator==(const EdgeCandidate& that) const
        {
            return tie(v0, v1) == tie(that.v0, that.v1);
        }
        bool operator<(const EdgeCandidate& that) const
        {
            return tie(v0, v1) < tie(that.v0, that.v1);
        }
    };
    vector<EdgeCandidate> edgeCandidates;
    void createEdgeCandidates();

    // Edge creation.
    void createEdges();
    void createEdgesMultithreaded();
    void createEdgesThreadFunction(uint64_t);

    // Prune removes all vertices that are not accessible from long
    // vertices in both directions.
    void prune();

    // Fill the Markov tables in the vertices.
    // These are needed to compute random paths.
    void fillMarkovTables();

    // Output.
    void write(const string& name) const;
    void writeCsv(const string& name) const;
    void writeVerticesCsv(const string& name) const;
    void writeEdgesCsv(const string& name) const;
    void writeGraphviz(const string& name) const;

    // Find the segment id of a vertex.
    uint64_t segmentId(vertex_descriptor v) const;

    // Random paths.
    // These functions find a random path starting at the given vertex.
    // Direction is 0 for forward and 1 backward.
    // The path stops when a long vertex is encountered.
    // Note these are paths in the ReadFollowing::Graph
    // but not in the AssemblyGraph.
    template<std::uniform_random_bit_generator RandomGenerator> void findRandomPath(
        vertex_descriptor, uint64_t direction,
        RandomGenerator&,
        Path&);

    // Python callable.
    void writeRandomPath(Segment, uint64_t direction);
    void writePathStatistics(Segment, uint64_t direction);
};
