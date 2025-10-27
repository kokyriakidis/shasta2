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

    class ReadFollowingVertex;
    class ReadFollowingEdge;

    using ReadFollowingBaseClass = boost::adjacency_list<
        boost::listS,
        boost::listS,
        boost::bidirectionalS,
        ReadFollowingVertex,
        ReadFollowingEdge>;

    namespace ReadFollowing {
        class Graph;
    }
}



class shasta::ReadFollowingVertex {
public:
    // A Segment is an edge of the AssemblyGraph.
    using Segment = AssemblyGraph::edge_descriptor;
    Segment segment;

    // The sequence length or estimated offset of this AssemblyGraph edge.
    uint64_t length = invalid<uint64_t>;

    ReadFollowingVertex(
        const AssemblyGraph&,
        AssemblyGraph::edge_descriptor);

};



class shasta::ReadFollowingEdge {
public:
    uint64_t coverage = 0;
    double jaccard = 0.;
};



class shasta::ReadFollowing::Graph : public ReadFollowingBaseClass {
public:
    Graph(const AssemblyGraph&);

private:
    const AssemblyGraph& assemblyGraph;

    // A Segment is an edge of the AssemblyGraph.
    using Segment = AssemblyGraph::edge_descriptor;

    // EXPOSE WHEN CODE STABILIZES.
    const uint64_t representativeRegionLength = 10;



    // Find appearances of OrientedReadIds in the initial/final
    // representative regions of each Segment.
    // Store them by OrientedReadId, indexed by OrientedReadId::getValue().
    // For the initial representative region, store the last appearance
    // (in journey order) of each OrientedReadId.
    // For the final representative region, store the first appearance
    // (in journey order) of each OrientedReadId.
    class AppearanceInfo {
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
    class Appearance : public AppearanceInfo {
    public:
        Segment segment;
        Appearance(const AppearanceInfo& appearanceInfo, Segment segment) :
            AppearanceInfo(appearanceInfo),
            segment(segment)
            {}
    };
    vector< vector<Appearance> > initialAppearances;
    vector< vector<Appearance> > finalAppearances;
    void findAppearances();



    // The number of initial/final Appearances in each Segment.
    std::map<Segment, uint64_t> initialAppearancesCount;
    std::map<Segment, uint64_t> finalAppearancesCount;
    void countAppearances();
    uint64_t getInitialAppearancesCount(Segment) const;
    uint64_t getFinalAppearancesCount(Segment) const;

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
