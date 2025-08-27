// Shasta.
#include "LocalAssemblyGraph.hpp"
#include "Anchor.hpp"
#include "graphvizToHtml.hpp"
#include "LocalSubgraph.hpp"
#include "tmpDirectory.hpp"
using namespace shasta;

// Boost libraries.
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>

// Standard library.
#include <fstream.hpp>



LocalAssemblyGraph::LocalAssemblyGraph(
    const AssemblyGraph& assemblyGraph,
    vector<AssemblyGraph::vertex_descriptor> startVertices,
    uint64_t maxDistance)
{
    *this = createLocalSubgraph<AssemblyGraph, LocalAssemblyGraph>(
        assemblyGraph, startVertices, true, true, maxDistance);
}



void LocalAssemblyGraph::writeHtml(
    ostream& html,
    const AssemblyGraph& assemblyGraph,
    uint64_t maxDistance) const
{
    // Write it out in Graphviz format.
    const string uuid = to_string(boost::uuids::random_generator()());
    const string dotFileName = tmpDirectory() + uuid + ".dot";
    writeGraphviz(dotFileName, assemblyGraph, maxDistance);


    // Display it in html in svg format.
    const double timeout = 30.;
    const string options = "-Nshape=rectangle -Gbgcolor=gray95";
    html << "<p>";
    graphvizToHtml(dotFileName, "dot", timeout, options, html);

}



void LocalAssemblyGraph::writeGraphviz(
    const string& fileName,
    const AssemblyGraph& assemblyGraph,
    uint64_t maxDistance) const
{
    ofstream dot(fileName);
    writeGraphviz(dot, assemblyGraph, maxDistance);
}



void LocalAssemblyGraph::writeGraphviz(
    ostream& dot,
    const AssemblyGraph& assemblyGraph,
    uint64_t maxDistance) const
{
    const LocalAssemblyGraph& localAssemblyGraph = *this;

    dot << "digraph LocalAssemblyGraph {\n";

    // Vertices.
    BGL_FORALL_VERTICES(lv, localAssemblyGraph, LocalAssemblyGraph) {
        const auto& localVertex = localAssemblyGraph[lv];
        const AssemblyGraph::vertex_descriptor v = localVertex.v;
        const uint64_t distance = localVertex.distance;
        const AssemblyGraphVertex& vertex = assemblyGraph[v];
        dot << vertex.id;

        // Begin vertex attributes.
        dot << "[";

        // Vertex label;
        dot << "label=\"" << vertex.id << "\\n" << anchorIdToString(vertex.anchorId);
        // dot << "\\n" << distance;
        dot << "\"";

        // Color by distance.
        if(distance == 0) {
            dot << " style=filled fillcolor=pink";
        } else if(distance == maxDistance) {
            dot << " style=filled fillcolor=cyan";
        } else {
            dot << " style=filled fillcolor=ivory";
        }

        // End vertex attributes.
        dot << "]";

        // End the line for this vertex.
        dot << ";\n";
    }



    // Edges.
    BGL_FORALL_EDGES(le, localAssemblyGraph, LocalAssemblyGraph) {
        const AssemblyGraph::edge_descriptor e = localAssemblyGraph[le].e;
        const AssemblyGraphEdge& edge = assemblyGraph[e];
        const uint64_t length = (edge.wasAssembled ? edge.sequenceLength() : edge.offset());
        const double penwidth = std::pow(double(length), 0.15);
        const vertex_descriptor lv0 = source(le, localAssemblyGraph);
        const vertex_descriptor lv1 = target(le, localAssemblyGraph);
        const AssemblyGraph::vertex_descriptor v0 = localAssemblyGraph[lv0].v;
        const AssemblyGraph::vertex_descriptor v1 = localAssemblyGraph[lv1].v;
        const AssemblyGraphVertex& vertex0 = assemblyGraph[v0];
        const AssemblyGraphVertex& vertex1 = assemblyGraph[v1];
        dot <<
            vertex0.id << "->" << vertex1.id <<
            "["
            "label=\"" << edge.id << "\\n" << length << "\""
            " penwidth=" << penwidth <<
            " color=navy"
            "];\n";
    }


    dot << "}\n";
}
