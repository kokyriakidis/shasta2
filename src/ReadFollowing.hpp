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
        class Appearance;

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



// Classes to describe an appearance of an OrientedReadIds in the initial/final
// representative regions of a Segment.
class shasta::ReadFollowing::Appearance {
public:
    Segment segment;
    uint64_t stepId;
    uint64_t offset; // To/from end/beginning of segment.

    OrientedReadId orientedReadId;
    uint32_t positionInJourney;
    uint32_t ordinal;
    uint32_t position;

    Appearance(
        Segment segment,
        uint64_t stepId,
        uint64_t offset,
        OrientedReadId orientedReadId,
        uint32_t positionInJourney,
        uint32_t ordinal,
        uint32_t position) :
        segment(segment),
        stepId(stepId),
        offset(offset),
        orientedReadId(orientedReadId),
        positionInJourney(positionInJourney),
        ordinal(ordinal),
        position(position)
        {}

    bool operator<(const Appearance& that) const
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

    vector<Appearance> initialAppearances;
    vector<Appearance> finalAppearances;

};



class shasta::ReadFollowing::Edge {
public:

    // Pairs of Appearances of the same OrientedReadId
    // in the final representative region of the source segment
    // and in the initial representative region of the target segment.
    vector< pair<Appearance, Appearance> > appearancePairs;

    uint64_t coverage() const
    {
        return appearancePairs.size();
    }

    double jaccard = 0.;
};



class shasta::ReadFollowing::Graph : public GraphBaseClass {
public:
    Graph(const AssemblyGraph&);

private:
    const AssemblyGraph& assemblyGraph;

    // Create vertices of the ReadFollowing graph.
    // Each vertex corresponds to a Segment.
    std::map<Segment, vertex_descriptor> vertexMap;
    void createVertices();

    // Create edges of the ReadFollowing graph.
    void createEdges();
    void writeEdgeDetails();

    // Enforce a minimum Jaccard when writing the graph.
    void writeGraph(double minJaccard) const;

    // Jaccard similarity for an edge.
    // Computed using the finalAppearancesCount of the source vertex
    // and the initialAppearancesCount of the target vertex.
    double jaccard(edge_descriptor) const;

public:
    void findPath(Segment, uint64_t direction, vector<vertex_descriptor>& path) const;
    void findForwardPath(Segment, vector<vertex_descriptor>& path) const;
    void findBackwardPath(Segment, vector<vertex_descriptor>& path) const;
    void writePath(Segment, uint64_t direction) const;
};
