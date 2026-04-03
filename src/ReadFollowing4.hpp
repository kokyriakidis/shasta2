#pragma once

/*****************************************************************

Read following in the AssemblyGraph.

*****************************************************************/

// Shasta.
#include "AssemblyGraph.hpp"
#include "MultithreadedObject.hpp"
#include "SegmentStepSupport.hpp"



namespace shasta2 {

    namespace ReadFollowing4 {

        class ReadFollower;

        class SearchGraph;
        class SearchGraphVertex;
        class SearchGraphEdge;
        using SearchGraphBaseClass = boost::adjacency_list<
            boost::listS,
            boost::listS,
            boost::bidirectionalS,
            SearchGraphVertex,
            SearchGraphEdge>;

        // A Segment is an edge of the AssemblyGraph.
        using Segment = AssemblyGraph::edge_descriptor;

        // EXPOSE WHEN CODE STABILIZES.
        const double logPThreshold = 0.;    // dB
        const double a = 3.;                // dB
        const double b = 10.;               // dB
    }

}



class shasta2::ReadFollowing4::SearchGraphVertex {
public:
    // A Segment is an edge of the AssemblyGraph.
    Segment segment;

    // The sequence length or estimated offset of this AssemblyGraph edge.
    uint64_t length = invalid<uint64_t>;

    // This is set for long vertices (length >= readFollowingSegmentLengthThreshold).
    bool isLong = false;

    SearchGraphVertex(Segment, uint64_t length, bool isLong);
};



class shasta2::ReadFollowing4::SearchGraphEdge {
public:
    uint64_t commonCount;
    uint64_t missingCount0;
    uint64_t missingCount1;

    SearchGraphEdge(
        uint64_t commonCount,
        uint64_t missingCount0,
        uint64_t missingCount1);

    // These are computed by the constructor from the three above.
    // logs are in dB.
    double logP;
    double logPForward;
    double logPBackward;
    double weight;

    double maxLogP() const;


    enum class Type {
        Bidirectional,
        Forward,
        Backward,
        Ambiguous
    };
    Type type() const;
};



class shasta2::ReadFollowing4::SearchGraph : public SearchGraphBaseClass {
public:

    // A map that gives the vertex_descriptor corresponding to each Segment.
    std::map<Segment, vertex_descriptor> vertexMap;

    // Create a vertex and update the vertexMap.
    void createVertex(Segment, uint64_t length, bool isLong);

    // Prune removes all vertices that are not accessible from long
    // vertices in both directions.
    void prune();

    // The vertex index map is needed to compute shortest paths.
    // It must be created when no more changes will be made to the graph.
    std::map<vertex_descriptor, uint64_t> vertexIndexMap;
    void createVertexIndexMap();

    void writeGraphviz(
        const AssemblyGraph&,
        const string& name,
        bool longGraph) const;

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


    void findAssemblyPaths(vector< vector<Segment> >& assemblyPaths) const;
    void writeAssemblyPaths(const vector< vector<Segment> >& assemblyPaths) const;
    void findAndWriteAssemblyPaths() const;

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

    // One graph for each read following direction (0=forward, 1=backward).
    // The graphs are identical except for edge directions.
    // This is necessary because I was not able to get boost::dijkstra_shortest_paths
    // to work in the reversed direction.
    array<SearchGraph, 2> graphs;

    // I also store a Graph which contains only vertices corresponding to long Segments.
    // This is not necessarily a subgraph of the above Graphs because it
    // allows "one-dorectional" edges in which missingCount0 is much less than missingCount1
    // or vice versa.
    SearchGraph longGraph;


    void createVertices();
    void createEdges();
    void createEdgesThreadFunction(uint64_t threadId);

    const AssemblyGraph& assemblyGraph;

private:
    bool isLong(Segment) const;
};
