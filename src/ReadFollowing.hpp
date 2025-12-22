#pragma once

// Read following in the AssemblyGraph.

// Shasta.
#include "AssemblyGraph.hpp"
#include "SegmentStepSupport.hpp"

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>

// Standard library.
#include <map>
#include <set>



namespace shasta2 {

    namespace ReadFollowing {
        class Graph;
        class Vertex;
        class Edge;

        using GraphBaseClass = boost::adjacency_list<
            boost::listS,
            boost::listS,
            boost::bidirectionalS,
            Vertex,
            Edge>;

        // A Segment is an edge of the AssemblyGraph.
        using Segment = AssemblyGraph::edge_descriptor;
    }
}



class shasta2::ReadFollowing::Vertex {
public:
    // A Segment is an edge of the AssemblyGraph.
    Segment segment;

    // The sequence length or estimated offset of this AssemblyGraph edge.
    uint64_t length = invalid<uint64_t>;

    Vertex(const AssemblyGraph&, Segment);

    vector<SegmentStepSupport> initialSupport;
    vector<SegmentStepSupport> finalSupport;

};



class shasta2::ReadFollowing::Edge {
public:
    Edge(const AssemblyGraph&, Segment, Segment);
    SegmentPairInformation segmentPairInformation;

    using Score = double;
    Score score;

    // For each edge v0->v1:
    // - isBest0 is set if this edge has the best score out of all out-edges of v0.
    // - isBest1 is set if this edge has the best score out of all in-edges of v1.
    // These are not maintained. They are set by setBestEdgeFlags().
    bool isBest0 = false;
    bool isBest1 = false;
};



class shasta2::ReadFollowing::Graph : public GraphBaseClass {
public:
    Graph(const AssemblyGraph&);

private:
    const AssemblyGraph& assemblyGraph;

    // Initial creation.
    std::map<Segment, vertex_descriptor> vertexMap;
    void createVertices();
    void createEdges();

    void computeEdgeScores();
    void setBestEdgeFlags();

    // Edge cleanup.
    void removeNegativeOffsetEdges();
    void removeLowCommonCountEdges(uint64_t minCommonCount);
    void removeLowCommonCorrectedJaccardEdges(double minCorrectedJaccard);

    // Prune sort leaves.
    void prune();
    bool pruneIteration();

    // Remove edges that don't have the best score at least one direction.
    void removeNonBestScoreEdges();

    // Output.
    void write(const string& name) const;
    void writeCsv(const string& name) const;
    void writeVerticesCsv(const string& name) const;
    void writeEdgesCsv(const string& name) const;
    void writeGraphviz(const string& name) const;

public:

    // Path following.
    // These functions find a best path starting at the given vertex and
    // ending if one of the forbiddenVertices is encountered.
    // Direction is 0 for forward and 1 backward.
    // Note these are paths in the ReadFollowing::Graph
    // but not in the AssemblyGraph.

    void findPath(
        vertex_descriptor, uint64_t direction,
        vector<vertex_descriptor>& path,
        const std::set<vertex_descriptor>& stopVertices) const;
    void findForwardPath(
        vertex_descriptor,
        vector<vertex_descriptor>& path,
        const std::set<vertex_descriptor>& stopVertices) const;
    void findBackwardPath(
        vertex_descriptor,
        vector<vertex_descriptor>& path,
        const std::set<vertex_descriptor>& stopVertices) const;

    // Find assembly paths.
    // These are best paths between vertices corresponding to long segments.
    void findPaths(vector< vector<Segment> >& assemblyPaths) const;

    // Python callable.
    void writePath(Segment, uint64_t direction) const;
    void writePaths() const;
};
