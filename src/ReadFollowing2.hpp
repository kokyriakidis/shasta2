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

        class Subgraph;
        class SubgraphVertex;
        class SubgraphEdge;
        using SubgraphBaseClass = boost::adjacency_list<
            boost::listS,
            boost::listS,
            boost::bidirectionalS,
            SubgraphVertex,
            SubgraphEdge>;

        using Path = vector<GraphBaseClass::vertex_descriptor>;

        class PathGraph;
        class PathGraphVertex;
        class PathGraphEdge;
        using PathGraphBaseClass = boost::adjacency_list<
            boost::listS,
            boost::listS,
            boost::bidirectionalS,
            PathGraphVertex,
            PathGraphEdge>;

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

        // Order by decreasing count, using segmentId to break ties.
        bool operator<(const RandomPathInfo& that) const
        {
            if(count > that.count) {
                return true;
            }
            if(count < that.count) {
                return false;
            }
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
    Graph(uint64_t iteration, const AssemblyGraph&);

    const AssemblyGraph& assemblyGraph;
    uint64_t iteration;

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

    // Find short vertices that, based on the stored random paths,
    // are reliably preceded/followed by a single vertex.
    // These will be used to fill in assembly paths.
    std::map< pair<vertex_descriptor, vertex_descriptor>, vector<vertex_descriptor> > randomPathsMap;
    void createRandomPathsMap();



    // Assembly paths.
    shared_ptr<PathGraph> createPathGraph();
    void findAssemblyPaths(vector< vector<Segment> >& assemblyPaths);
    void writeAssemblyPaths(const vector< vector<Segment> >& assemblyPaths) const;

    // Find an assembly path between two long vertices.
    // This uses the randomPathsMap to locate usable short vertices.
    // The Path does not include the segments corresponding to v0 and v1.
    // This can fail, in which case it returns false and an empty assembly path.
    bool findAssemblyPath(
        vertex_descriptor v0,
        vertex_descriptor v1,
        vector<Segment>&) const;



    // Python callable.
    void writeRandomPath(Segment, uint64_t direction);
    void writePathStatistics(Segment, uint64_t direction);
    void findAndWriteAssemblyPaths();

    class OrderById {
    public:
        bool operator() (const vertex_descriptor v0, const vertex_descriptor v1) const
        {
            return graph.segmentId(v0) < graph.segmentId(v1);
        }
        OrderById(const Graph& graph) : graph(graph) {}
        const Graph& graph;

    };
    const OrderById orderById;
};



class shasta2::ReadFollowing2::SubgraphVertex {
public:
    Segment segment;

    // Used by approximateTopologicalSort.
    uint64_t rank = 0;
    uint64_t color = 0;
};



class shasta2::ReadFollowing2::SubgraphEdge {
public:
    double correctedJaccard = 0.;
    bool isDagEdge = false;
};



// A subgraph of the Graph.
// It only stores a subset of the information stored in the Graph.
class shasta2::ReadFollowing2::Subgraph : public SubgraphBaseClass {
public:
    Subgraph(
        const Graph&,
        const std::set<Graph::vertex_descriptor, Graph::OrderById>& vertices);
    const Graph& graph;
    uint64_t segmentId(vertex_descriptor) const;
    void makeAcyclic();
    void writeGraphviz(const string& name) const;
};



// Each PathGraph vertex corresponds to a long segment.
class shasta2::ReadFollowing2::PathGraphVertex {
public:
    Segment segment;
};



class shasta2::ReadFollowing2::PathGraphEdge {
public:

    // If segment0 and segment1 are the source and target segments
    // for this PathGraph edge:
    // - randomPathCount[0] contains the number of forward paths starting from segment0
    //   that ended at segment1.
    // - randomPathCount[1] contains the number of backward paths starting from segment1
    //   that ended at segment0.
    array<uint64_t, 2> randomPathCount = {0, 0};

    // A possible assembly path between the long segments corresponding
    // to the vertices of this PathGraphEdge.
    vector<Segment> assemblyPath;
};



// Class used to store paths between long segments.
class shasta2::ReadFollowing2::PathGraph : public PathGraphBaseClass {
public:

    // This creates an empty PathGraph.
    PathGraph(const Graph&);

    // This uses the Graph to create vertices aned edges.
    void create();

    void findAssemblyPathsOnEdges();

    // This returns all the non-trivial connected components
    // of the PathGraph (that is, the ones with at least two vertices).
    vector< shared_ptr<PathGraph> > findConnectedComponents();

    // This computes an assembly path, assuming it is working
    // on a PathGraph with a single connected component.
    void findAssemblyPath(vector<Segment>&);

    uint64_t segmentId(vertex_descriptor) const;
private:
    const Graph& graph;
    void createVertices();
    void createEdges();

    // Find assembly paths on all edges, then remove the
    // edges for which an assembly path could not be found.

    // Map segments to vertices.
    std::map<Segment, vertex_descriptor> vertexMap;

public:
    // Graphviz output.
    void writeGraphviz(const string& name) const;
    void writeGraphviz(ostream&) const;
};
