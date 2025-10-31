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



namespace shasta {

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



class shasta::ReadFollowing::Vertex {
public:
    // A Segment is an edge of the AssemblyGraph.
    Segment segment;

    // The sequence length or estimated offset of this AssemblyGraph edge.
    uint64_t length = invalid<uint64_t>;

    Vertex(const AssemblyGraph&, Segment);

    vector<SegmentStepSupport> initialSupport;
    vector<SegmentStepSupport> finalSupport;

};



class shasta::ReadFollowing::Edge {
public:
    Edge(const AssemblyGraph&, Segment, Segment);
    SegmentPairInformation segmentPairInformation;

    // For each edge v0->v1:
    // - isLowestOffset0 is set if this edge has the lowest offset out of all out-edges of v0.
    // - isLowestOffset1 is set if this edge has the lowest offset out of all in-edges of v1.
    // These are not maintained. They are set by setLowestOffsetFlags().
    bool isLowestOffset0 = false;
    bool isLowestOffset1 = false;
};



class shasta::ReadFollowing::Graph : public GraphBaseClass {
public:
    Graph(const AssemblyGraph&);

private:
    const AssemblyGraph& assemblyGraph;

    // Initial creation.
    std::map<Segment, vertex_descriptor> vertexMap;
    void createVertices();
    void createEdges();

    // Remove edges with negative offset.
    void removeNegativeOffsetEdges();

    // Remove weak edges.
    void removeLowCommonCountEdges(uint64_t minCommonCount);
    void removeLowCommonCorrectedJaccardEdges(double minCorrectedJaccard);

    // Prune sort leaves.
    void prune(uint64_t minimumLength);
    bool pruneIteration(uint64_t minimumLength);

    // Remove edges that have both isLowestOffset0 and isLowestOffset1 set to false.
    void removeNonLowestOffsetEdges();

    void setLowestOffsetFlags();

    void write(const string& name);
    void writeCsv(const string& name) const;
    void writeVerticesCsv(const string& name) const;
    void writeEdgesCsv(const string& name) const;
    void writeGraphviz(const string& name);

public:

    // Find a minimum offset path starting at the given vertex and
    // ending if one of the terminalVertices is encountered.
    // Direction is 0 for forward and 1 backward.
    void findPath(
        Segment, uint64_t direction,
        vector<vertex_descriptor>& path,
        const std::set<vertex_descriptor>& terminalVertices) const;
    void findForwardPath(
        Segment,
        vector<vertex_descriptor>& path,
        const std::set<vertex_descriptor>& terminalVertices) const;
    void findBackwardPath(
        Segment,
        vector<vertex_descriptor>& path,
        const std::set<vertex_descriptor>& terminalVertices) const;

    // Python callable.
    void writePath(Segment, uint64_t direction) const;
};
