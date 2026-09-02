#pragma once

// A SegmentGraph is a directed graph in which each vertex
// corresponds to an AssemblyGraph edge, that is, a Segment.
// The meaning of a SegmentGraph edge is dependent on what the
// SegmentGraph is used for.
// One can construct a SegmentGraph which is exactly a line graph
// of the AssemblyGraph (https://en.wikipedia.org/wiki/Line_graph).
// In this case an edge s0->s1
// means that, in the AssemblyGraph, the target vertex of s0
// coincides with the source vertex of s1.

#include "AssemblyGraphBaseClass.hpp"
#include "SegmentStepSupport.hpp"

namespace shasta2 {

    class SegmentGraph;

    using SegmentGraphBaseClass = boost::adjacency_list<
        boost::listS,
        boost::listS,
        boost::bidirectionalS,
        Segment,
        SegmentPairInformation>;

    class AssemblyGraph;
}



class shasta2::SegmentGraph: public SegmentGraphBaseClass {
public:

    // Construct a SegmentGraph which is the line graph of an AssemblyGraph.
    // (https://en.wikipedia.org/wiki/Line_graph).
    // But if createEdges is false, only the vertices are created.
    SegmentGraph(const AssemblyGraph&, bool createEdges);

    // The vertexMap gives the vertex descriptor corresponding to a Segment.
    std::map<Segment, vertex_descriptor> vertexMap;

    const AssemblyGraph& assemblyGraph;

    void addEdge(Segment, Segment);
    void addEdge(Segment, Segment, const SegmentPairInformation&);

    void writeGraphviz(const string& fileName) const;
    void writeGraphviz(ostream&) const;
};
