#include "SegmentGraph.hpp"
#include "AssemblyGraph.hpp"
using namespace shasta2;

#include "fstream.hpp"



// Construct a SegmentGraph which is the line graph of an AssemblyGraph.
// (https://en.wikipedia.org/wiki/Line_graph).
// But if createEdges is false, only the vertices are created.
SegmentGraph::SegmentGraph(
    const AssemblyGraph& assemblyGraph, bool createEdges) :
    assemblyGraph(assemblyGraph)
{
    SegmentGraph& segmentGraph = *this;

    // Create the vertices.
    BGL_FORALL_EDGES(segment, assemblyGraph, AssemblyGraph) {
        const vertex_descriptor v = add_vertex(segment, segmentGraph);
        vertexMap.insert({segment, v});
    }

    // Create the edges, if requested.
    if(not createEdges) {
        return;
    }
    BGL_FORALL_VERTICES(av, assemblyGraph, AssemblyGraph) {
        BGL_FORALL_INEDGES(av, segment0, assemblyGraph, AssemblyGraph) {
            BGL_FORALL_OUTEDGES(av, segment1, assemblyGraph, AssemblyGraph) {
               addEdge(segment0, segment1);
            }
        }
    }
}



void SegmentGraph::addEdge(
    Segment s0,
    Segment s1)
{
    add_edge(vertexMap.at(s0), vertexMap.at(s1), *this);
}



void SegmentGraph::addEdge(
    Segment s0,
    Segment s1,
    const SegmentPairInformation& segmentPairInformation)
{
    add_edge(vertexMap.at(s0), vertexMap.at(s1), segmentPairInformation, *this);
}



void SegmentGraph::writeGraphviz(const string& fileName) const
{
    ofstream dot(fileName);
    writeGraphviz(dot);
}



void SegmentGraph::writeGraphviz(ostream& dot) const
{
    const SegmentGraph segmentGraph = *this;
    dot << "digraph SegmentGraph {\n";

    // Only write the edges. this way the isolated Segments are not includes.
    BGL_FORALL_EDGES(e, segmentGraph, SegmentGraph) {
        const SegmentPairInformation segmentPairInformation = segmentGraph[e];
        const vertex_descriptor v0 = source(e, segmentGraph);
        const vertex_descriptor v1 = target(e, segmentGraph);
        const Segment s0 = segmentGraph[v0];
        const Segment s1 = segmentGraph[v1];
        dot << assemblyGraph.id(s0) << "->" <<
            assemblyGraph.id(s1) <<
            "[penwidth=\"" << 0.2 * double(segmentPairInformation.commonCount) << "\""
            " tooltip=\"" << segmentPairInformation.commonCount << "/" <<
            segmentPairInformation.missing() << "/" <<
            segmentPairInformation.segmentOffset <<
            "\""
            "]"
            ";\n";

    }

    dot << "}\n";
}
