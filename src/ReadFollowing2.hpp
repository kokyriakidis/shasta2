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

    namespace ReadFollowing2 {

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



class shasta2::ReadFollowing2::Vertex {
public:
    // A Segment is an edge of the AssemblyGraph.
    Segment segment;

    // The sequence length or estimated offset of this AssemblyGraph edge.
    uint64_t length = invalid<uint64_t>;

    // This is set for long vertices (length >= readFollowingSegmentLengthThreshold).
    bool isLongVertex = false;

    double coverage = 0.;

    Vertex(const AssemblyGraph&, Segment);

    vector<SegmentStepSupport> initialSupport;
    vector<SegmentStepSupport> finalSupport;

    // Vectors used for the generation of random paths.
    // Each child/parent vertex of this vertex appears
    // a number of times equal to its commonCount.
    // This way, the random paths are biased by commonCount.
    // These are filled in by fillConnectivity.
    vector<GraphBaseClass::vertex_descriptor> children;
    vector<GraphBaseClass::vertex_descriptor> parents;

    // For each starting vertex we generate pathCount random paths
    // in each direction. We store the terminal vertices reached by
    // each of the random paths.
    // For a forward path (direction = 0), the terminal vertex is the last vertex.
    // For a backward path (direction = 0), the terminal vertex is the first vertex.
    class RandomPathInfo {
    public:
        GraphBaseClass::vertex_descriptor v;
        uint64_t count;
        uint64_t segmentId;
        bool operator<(const RandomPathInfo& that) const
        {
            return segmentId < that.segmentId;
        }
    };
    array<vector<RandomPathInfo>, 2> randomPathInfos;


};



class shasta2::ReadFollowing2::Edge {
public:
    Edge(const AssemblyGraph&, Segment, Segment);
    SegmentPairInformation segmentPairInformation;
};



class shasta2::ReadFollowing2::Graph :
    public GraphBaseClass,
    public MultithreadedObject<Graph> {
public:
    Graph(const AssemblyGraph&);

    // EXPOSE WHEN CODE STABILIZES.
    static const uint64_t pathCount = 100;
    static constexpr double minPathCountFraction = 0.4;

    const AssemblyGraph& assemblyGraph;

    // Vertex creation.
    std::map<Segment, vertex_descriptor> vertexMap;
    void createVertices();
    std::set<vertex_descriptor> longVertices;

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

    // Prune removes all vertices that are not accessible from "long"
    // vertices in both directions.
    void prune();

    // This fills in the outEdges and inEdges vectors
    // of all vertices, which are then used to generate
    // random paths.
    void fillConnectivity();

    // Output.
    void write(const string& name) const;
    void writeCsv(const string& name) const;
    void writeVerticesCsv(const string& name) const;
    void writeEdgesCsv(const string& name) const;
    void writeGraphviz(const string& name) const;

    uint64_t segmentId(vertex_descriptor v) const;



    // Random paths.
    // These functions find a random path starting at the given vertex.
    // Direction is 0 for forward and 1 backward.
    // The path stops when a long vertex is encountered.
    // Note these are paths in the ReadFollowing1::Graph
    // but not in the AssemblyGraph.
    template<std::uniform_random_bit_generator RandomGenerator> void findRandomPath(
        vertex_descriptor, uint64_t direction,
        RandomGenerator&,
        Path&);
    template<std::uniform_random_bit_generator RandomGenerator> void findRandomForwardPath(
        vertex_descriptor,
        RandomGenerator&,
        Path&);
    template<std::uniform_random_bit_generator RandomGenerator> void findRandomBackwardPath(
        vertex_descriptor,
        RandomGenerator&,
        Path&);

    // Create pathCount paths in each starting direction and for each starting vertex.
    // Store in each vertex the terminal vertices of the random paths.
    void findRandomPaths();

    // Python callable.
    void writeRandomPath(Segment, uint64_t direction);
    void writePathStatistics(Segment, uint64_t direction);
};

