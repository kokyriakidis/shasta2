// Shasta.
#include "LocalAssemblyGraph.hpp"
#include "Anchor.hpp"
#include "approximateTopologicalSort.hpp"
#include "graphvizToHtml.hpp"
#include "LocalSubgraph.hpp"
#include "orderPairs.hpp"
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
    uint64_t maxDistance)
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
    uint64_t maxDistance)
{
    ofstream dot(fileName);
    writeGraphviz(dot, assemblyGraph, maxDistance);
}



void LocalAssemblyGraph::writeGraphviz(
    ostream& dot,
    const AssemblyGraph& assemblyGraph,
    uint64_t maxDistance)
{
    LocalAssemblyGraph& localAssemblyGraph = *this;

    // Compute an approximate topological sort.
    vector< pair<edge_descriptor, double> > edgesWithCoverage;
    BGL_FORALL_EDGES(le, localAssemblyGraph, LocalAssemblyGraph) {
        const AssemblyGraph::edge_descriptor e = localAssemblyGraph[le].e;
        edgesWithCoverage.emplace_back(le, assemblyGraph[e].averageCoverage());
    }
    std::ranges::sort(edgesWithCoverage, OrderPairsBySecondOnlyGreater<edge_descriptor, double>());
    vector<edge_descriptor> edgesSortedByCoverage;
    for(const auto& p: edgesWithCoverage) {
        edgesSortedByCoverage.push_back(p.first);
    }
    approximateTopologicalSort(localAssemblyGraph, edgesSortedByCoverage);

    // Gather vertices sorted by rank.
    vector< pair<vertex_descriptor, uint64_t> > verticesWithRank;
    BGL_FORALL_VERTICES(lv, localAssemblyGraph, LocalAssemblyGraph) {
        verticesWithRank.emplace_back(lv, localAssemblyGraph[lv].rank);
    }
    std::ranges::sort(verticesWithRank, OrderPairsBySecondOnly<vertex_descriptor, uint64_t>());



    dot << "digraph LocalAssemblyGraph {\n";

    // Vertices. Written in rank order for better display.
    for(const auto& p: verticesWithRank) {
        const vertex_descriptor lv = p.first;
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
        const LocalAssemblyGraphEdge localEdge = localAssemblyGraph[le];
        const AssemblyGraph::edge_descriptor e = localEdge.e;
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
            " color=" << (localEdge.isDagEdge ? "navy" : "red") <<
            "];\n";
    }


    dot << "}\n";
}
