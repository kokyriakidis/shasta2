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
        return GraphBaseClass::null_vertex();
    }



    // For each starting vertex we generate pathCount random paths
    // in each direction. We store the terminal vertices reached by
    // each of the random paths.
    // For a forward path (direction = 0), the terminal vertex is the last vertex.
    // For a backward path (direction = 1), the terminal vertex is the first vertex.
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



class shasta2::ReadFollowing3::Edge {
public:
    Edge(const AssemblyGraph&, Segment, Segment);
    SegmentPairInformation segmentPairInformation;

    // Unnormalized probability of this edge, for random paths.
    double pUnnormalized = 0.;

    // The weight is the inverse of pUnnormalized.
    double weight;
};



class shasta2::ReadFollowing3::Graph :
    public GraphBaseClass,
    public MultithreadedObject<Graph> {
public:
    Graph(const AssemblyGraph&, bool createEmpty = false);

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

    // Create pathCount paths in each starting direction and for each starting vertex.
    // Store in each vertex the terminal vertices of the random paths.
    void findRandomPaths();

    // Find short vertices that, based on the stored random paths,
    // are reliably preceded/followed by a single vertex.
    // These will be used to fill in assembly paths.
    std::map< pair<vertex_descriptor, vertex_descriptor>, vector<vertex_descriptor> > randomPathsMap;
    void createRandomPathsMap();

    // Assembly paths.
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



    // This finds a shortest path starting at v0 and ending at a long vertex,
    // with path length defined by Edge::weight = 1/Edge::pUnnormalized.
    // So the shortest path prefers edges with high pUnnormalized.
    void findShortestPath(
        vertex_descriptor v0,   // The start vertex.
        uint64_t direction,     // 0 = forward, 1 = backward
        vector<vertex_descriptor>& path
        ) const;
    void findShortestPathForward(
        vertex_descriptor v0,   // The start vertex.
        vector<vertex_descriptor>& path
        ) const;
    void findShortestPathBackward(
        vertex_descriptor v0,   // The start vertex.
        vector<vertex_descriptor>& path
        ) const;
    void findAndWriteShortestPath(Segment, uint64_t direction) const; // Python callable

    // The vertex index map is needed to compute shortest paths.
    // It must be created when no more changes will be made to the graph.
    std::map<vertex_descriptor, uint64_t> vertexIndexMap;
    void createVertexIndexMap();


    // For findShortestPathBackward we also need a reversed version of the Graph.
    // I tried using boost::reverse_graph but it was problematic to use with
    // dijkstra_shortest_paths.
    shared_ptr<Graph> reversedGraphPointer;
    void createReversedGraph();
};



class shasta2::ReadFollowing3::SubgraphVertex {
public:
    Segment segment;

    // Used by approximateTopologicalSort.
    uint64_t rank = 0;
    uint64_t color = 0;
};



class shasta2::ReadFollowing3::SubgraphEdge {
public:
    double correctedJaccard = 0.;
    bool isDagEdge = false;
};



// A subgraph of the Graph.
// It only stores a subset of the information stored in the Graph.
class shasta2::ReadFollowing3::Subgraph : public SubgraphBaseClass {
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
class shasta2::ReadFollowing3::PathGraphVertex {
public:
    Segment segment;
};



class shasta2::ReadFollowing3::PathGraphEdge {
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

    // The shortest paths found in each direction.
    // At least one of them is non-empty.
    array< vector<Segment>, 2> assemblyPaths;
    bool isBidirectional() const
    {
        return (not assemblyPaths[0].empty()) and (not assemblyPaths[1].empty());
    }

    const vector<Segment>& assemblyPathToUse() const
    {
        if(not assemblyPaths[0].empty()) {
            return assemblyPaths[0];
        }
        if(not assemblyPaths[1].empty()) {
            return assemblyPaths[1];
        }
        SHASTA2_ASSERT(0);
    }
};



// Class used to store paths between long segments.
class shasta2::ReadFollowing3::PathGraph : public PathGraphBaseClass {
public:

    // This creates an empty PathGraph.
    PathGraph(const Graph&);

    // This uses the Graph to create vertices aned edges.
    void create();

    void removeWeakEdges();

    // Find assembly paths on all edges, then remove the
    // edges for which an assembly path could not be found.
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
