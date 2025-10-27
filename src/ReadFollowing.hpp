#pragma once

// Read following in the AssemblyGraph.

// Shasta.
#include "AssemblyGraph.hpp"

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
        class AppearanceInfo;
        class Appearance;

        using GraphBaseClass = boost::adjacency_list<
            boost::listS,
            boost::listS,
            boost::bidirectionalS,
            Vertex,
            Edge>;
        using Segment = AssemblyGraph::edge_descriptor;
    }
}



// Classes to describe an appearance of an OrientedReadIds in the initial/final
// representative regions of a Segment.
class shasta::ReadFollowing::AppearanceInfo {
public:
    uint32_t positionInJourney;
    uint32_t ordinal;
    uint32_t position;
    uint64_t stepId;

    // For an AppearanceInfo in the initial representative region,
    // store the base offset between this appearance and the end of the segment.
    // For an AppearanceInfo in the final representative region,
    // store the base offset between the beginning of the segment and this appearance.
    uint64_t offset;

    AppearanceInfo(
        uint32_t positionInJourney,
        uint32_t ordinal,
        uint32_t position,
        uint64_t stepId,
        uint64_t offset) :
        positionInJourney(positionInJourney),
        ordinal(ordinal),
        position(position),
        stepId(stepId),
        offset(offset)
    {}
    bool operator<(const AppearanceInfo& that) const
    {
        return positionInJourney < that.positionInJourney;
    }
};



class shasta::ReadFollowing::Vertex {
public:
    // A Segment is an edge of the AssemblyGraph.
    Segment segment;

    // The sequence length or estimated offset of this AssemblyGraph edge.
    uint64_t length = invalid<uint64_t>;

    Vertex(
        const AssemblyGraph&,
        AssemblyGraph::edge_descriptor);

    vector<AppearanceInfo> initialAppearances;
    vector<AppearanceInfo> finalAppearances;

};



class shasta::ReadFollowing::Edge {
public:
    uint64_t coverage = 0;
    double jaccard = 0.;
};



class shasta::ReadFollowing::Appearance : public AppearanceInfo {
public:
    Segment segment;
    Appearance(const AppearanceInfo& appearanceInfo, Segment segment) :
        AppearanceInfo(appearanceInfo),
        segment(segment)
        {}
};



class shasta::ReadFollowing::Graph : public GraphBaseClass {
public:
    Graph(const AssemblyGraph&);

private:
    const AssemblyGraph& assemblyGraph;

    // EXPOSE WHEN CODE STABILIZES.
    const uint64_t representativeRegionLength = 10;

    // Appearance of OrientedReadIds in the initial/final
    // representative regions of segments.
    // Indexed by OrientedReadId::getValue().
    vector< vector<Appearance> > initialAppearances;
    vector< vector<Appearance> > finalAppearances;
    void findAppearances();

    // Create vertices of the ReadFollowing graph.
    // Each vertex corresponds to a Segment, but not
    // all Segments generate a vertex.
    std::map<Segment, vertex_descriptor> vertexMap;
    void createVertices();

    // Create edges of the ReadFollowing graph.
    void createEdges();

    // Enforce a minimum Jaccard when writing the graph.
    void writeGraph(double minJaccard) const;

    // Jaccard similarity for an EdgePairsGraph edge.
    // Computed using the finalAppearancesCount of the source vertex
    // and the initialAppearancesCount of the target vertex.
    double jaccard(edge_descriptor) const;

public:
    void findPath(Segment, uint64_t direction, vector<vertex_descriptor>& path) const;
    void findForwardPath(Segment, vector<vertex_descriptor>& path) const;
    void findBackwardPath(Segment, vector<vertex_descriptor>& path) const;
    void writePath(Segment, uint64_t direction) const;
};
