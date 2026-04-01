#pragma once

/*****************************************************************

Read following in the AssemblyGraph.

This creates two directed graphs, one for each read following
direction (0=forward, 1=backward), in which each vertex
represents an edge of the AssemblyGraph, here called a
Segment using GFA terminology.

*****************************************************************/

// Shasta.
#include "AssemblyGraph.hpp"
#include "MultithreadedObject.hpp"
#include "SegmentStepSupport.hpp"



namespace shasta2 {

    namespace ReadFollowing4 {

        class Graph;
        class Vertex;
        class Edge;
        using GraphBaseClass = boost::adjacency_list<
            boost::listS,
            boost::listS,
            boost::bidirectionalS,
            Vertex,
            Edge>;

        class ReadFollower;


        // A Segment is an edge of the AssemblyGraph.
        using Segment = AssemblyGraph::edge_descriptor;

        // EXPOSE WHEN CODE STABILIZES.
        const double logPThreshold = 0.;    // dB
        const double a = 3.;                // dB
        const double b = 10.;               // dB
    }

}



class shasta2::ReadFollowing4::Vertex {
public:
    // A Segment is an edge of the AssemblyGraph.
    Segment segment;

    // The sequence length or estimated offset of this AssemblyGraph edge.
    uint64_t length = invalid<uint64_t>;

    // This is set for long vertices (length >= readFollowingSegmentLengthThreshold).
    bool isLong = false;

    Vertex(const AssemblyGraph&, Segment);
    Vertex() {}
};



// The Edge is assigned a probability based on missing/common
// counts for the SegmentPairInformation for the two vertices.
// The probability then is used as a weight for shortest paths.
class shasta2::ReadFollowing4::Edge {
public:
    uint64_t commonCount;
    uint64_t missingCount;
    double logP; // In decibels (dB)
    double weight;
};



class shasta2::ReadFollowing4::Graph : public GraphBaseClass {
public:

    // A map that gives the vertex_descriptor corresponding to each Segment.
    std::map<Segment, vertex_descriptor> vertexMap;

    // Create a vertex and update the vertexMap.
    void createVertex(const AssemblyGraph&, Segment);

    // Prune removes all vertices that are not accessible from long
    // vertices in both directions.
    void prune();

    // The vertex index map is needed to compute shortest paths.
    // It must be created when no more changes will be made to the graph.
    std::map<vertex_descriptor, uint64_t> vertexIndexMap;
    void createVertexIndexMap();

    void writeGraphviz(const AssemblyGraph&, const string& name) const;

    void findShortestPath(Segment, vector<Segment>&) const;
};



class shasta2::ReadFollowing4::ReadFollower :
    public MultithreadedObject<ReadFollower> {

public:
    ReadFollower(const AssemblyGraph&);

    // This finds a shortest path starting at segment0 and ending at a long Segment,
    // with path length defined by Edge::weight.
    void findShortestPath(
        Segment segment0,
        uint64_t direction,     // 0 = forward, 1 = backward
        vector<Segment>& path
        ) const;
    void findShortestPathForward(
        Segment segment0,
        vector<Segment>& path
        ) const;
    void findShortestPathBackward(
        Segment segment0,
        vector<Segment>& path
        ) const;
    void findAndWriteShortestPath(Segment, uint64_t direction) const; // Python callable

private:
    const AssemblyGraph& assemblyGraph;

    // Initial and final support for each Segment.
    std::map<Segment, vector<SegmentStepSupport> > initialSupportMap;
    std::map<Segment, vector<SegmentStepSupport> > finalSupportMap;
    void fillSupportMaps();

    // Segment pairs (segment0, segment1) such that the final support
    // of segment0 shares at least readFollowingMinCommonCount
    // OrientedReadIds with the initial support of segment1.
    // Each of them is a candidate edge for the two Graphs.
    vector< pair<Segment, Segment> > segmentPairs;
    void findSegmentPairs();

    // One graph for each read following direction
    // (0=forward, 1=backward).
    array<Graph, 2> graphs;

    void createVertices();
    void createEdges();
    void createEdgesThreadFunction(uint64_t threadId);

};

