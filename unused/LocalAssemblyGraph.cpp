// Shasta.
#include "LocalAssemblyGraph.hpp"
#include "Anchor.hpp"
#include "approximateTopologicalSort.hpp"
#include "deduplicate.hpp"
#include "graphvizToHtml.hpp"
#include "LocalSubgraph.hpp"
#include "orderPairs.hpp"
#include "shastaLapack.hpp"
#include "tmpDirectory.hpp"
using namespace shasta2;

// Boost libraries.
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>

// Standard library.
#include "fstream.hpp"



// Constructor from a set of AssemblyGraph start vertices and maximum distance.
LocalAssemblyGraph::LocalAssemblyGraph(
    const AssemblyGraph& assemblyGraph,
    vector<AssemblyGraph::vertex_descriptor> startVertices,
    uint64_t maxDistance)
{
    *this = createLocalSubgraph<AssemblyGraph, LocalAssemblyGraph>(
        assemblyGraph, startVertices, true, true, maxDistance);
    approximateTopologicalSort(assemblyGraph);
}



// Constructor from a set of AssemblyGraph edges.
LocalAssemblyGraph::LocalAssemblyGraph(
    const AssemblyGraph& assemblyGraph,
    vector<AssemblyGraph::edge_descriptor> edges)
{
    LocalAssemblyGraph& localAssemblyGraph = *this;

    // Create the vertices.
    std::map<AssemblyGraph::vertex_descriptor, vertex_descriptor> vertexMap;
    for(const AssemblyGraph::edge_descriptor e: edges) {

        const AssemblyGraph::vertex_descriptor v0 = source(e, assemblyGraph);
        if(not vertexMap.contains(v0)) {
            vertexMap.insert(make_pair(v0, add_vertex(LocalAssemblyGraphVertex(v0), localAssemblyGraph)));
        }

        const AssemblyGraph::vertex_descriptor v1 = target(e, assemblyGraph);
        if(not vertexMap.contains(v1)) {
            vertexMap.insert(make_pair(v1, add_vertex(LocalAssemblyGraphVertex(v1), localAssemblyGraph)));
        }
    }



    // Create the edges.
    for(const AssemblyGraph::edge_descriptor e: edges) {

        const AssemblyGraph::vertex_descriptor v0 = source(e, assemblyGraph);
        const auto it0 = vertexMap.find(v0);
        SHASTA2_ASSERT(it0 != vertexMap.end());
        const vertex_descriptor lv0 = it0->second;

        const AssemblyGraph::vertex_descriptor v1 = target(e, assemblyGraph);
        const auto it1 = vertexMap.find(v1);
        SHASTA2_ASSERT(it1 != vertexMap.end());
        const vertex_descriptor lv1 = it1->second;

        add_edge(lv0, lv1, LocalAssemblyGraphEdge(e), localAssemblyGraph);
    }

    approximateTopologicalSort(assemblyGraph);
}



void LocalAssemblyGraph::approximateTopologicalSort(const AssemblyGraph& assemblyGraph)
{
    LocalAssemblyGraph& localAssemblyGraph = *this;

    // Gather edges sorted by decreasing coverage.
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

    // Do the approximate topological sort.
    shasta2::approximateTopologicalSort(localAssemblyGraph, edgesSortedByCoverage);

    // Gather vertices sorted by rank.
    vector< pair<vertex_descriptor, uint64_t> > verticesWithRank;
    BGL_FORALL_VERTICES(lv, localAssemblyGraph, LocalAssemblyGraph) {
        verticesWithRank.emplace_back(lv, localAssemblyGraph[lv].rank);
    }
    std::ranges::sort(verticesWithRank, OrderPairsBySecondOnly<vertex_descriptor, uint64_t>());

    verticesByRank.clear();
    for(const auto& p: verticesWithRank) {
        verticesByRank.push_back(p.first);
    }

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

    // Vertices. Written in rank order for better display.
    for(const vertex_descriptor lv: verticesByRank) {
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



void LocalAssemblyGraph::analyze(
    ostream& html,
    const AssemblyGraph& assemblyGraph) const
{
    const LocalAssemblyGraph localAssemblyGraph = *this;

    // Gather the edges and sort them by id.
    vector<AssemblyGraph::edge_descriptor> edgesSortedById;
    BGL_FORALL_EDGES(le, localAssemblyGraph, LocalAssemblyGraph) {
        const AssemblyGraph::edge_descriptor e = localAssemblyGraph[le].e;
        edgesSortedById.push_back(e);
    }
    std::ranges::sort(edgesSortedById, assemblyGraph.orderById);
    const uint64_t edgeCount = edgesSortedById.size();

    // Gather the OrientedReadids present in all the edges.
    vector<OrientedReadId> orientedReadIds;
    for(const AssemblyGraph::edge_descriptor e: edgesSortedById) {
        const AssemblyGraphEdge& edge = assemblyGraph[e];
        for(const AssemblyGraphEdgeStep& step: edge) {
            const AnchorPair anchorPair = step.anchorPair;
            std::ranges::copy(anchorPair.orientedReadIds, back_inserter(orientedReadIds));
        }
    }
    deduplicate(orientedReadIds);
    const uint64_t orientedReadCount = orientedReadIds.size();



    // Count how many times each OrientedReadId appears in each edge.
    using IntegerMatrix = boost::numeric::ublas::matrix<uint64_t, boost::numeric::ublas::column_major>;
    IntegerMatrix integerMatrix(orientedReadCount, edgeCount, 0);
    for(uint64_t iEdge=0; iEdge<edgeCount; iEdge++) {
        const AssemblyGraph::edge_descriptor e = edgesSortedById[iEdge];
        const AssemblyGraphEdge& edge = assemblyGraph[e];
        for(const AssemblyGraphEdgeStep& step: edge) {
            const AnchorPair anchorPair = step.anchorPair;
            for(const OrientedReadId orientedReadId: anchorPair.orientedReadIds) {
                const auto it = std::ranges::lower_bound(orientedReadIds, orientedReadId);
                SHASTA2_ASSERT(it != orientedReadIds.end());
                SHASTA2_ASSERT(*it == orientedReadId);
                const uint64_t iOrientedRead = it - orientedReadIds.begin();
                ++integerMatrix(iOrientedRead, iEdge);
            }
        }
    }


    if(html) {
        html <<
            "<h3>Occurrences of oriented reads in edges of this local assembly graph</h3>"
            "<table><tr><td>";
        for(uint64_t iEdge=0; iEdge<edgeCount; iEdge++) {
            html << "<th>" << assemblyGraph[edgesSortedById[iEdge]].id;
        }
        for(uint64_t iOrientedRead=0; iOrientedRead<orientedReadCount; iOrientedRead++) {
            html << "<tr><th>" << orientedReadIds[iOrientedRead];
            for(uint64_t iEdge=0; iEdge<edgeCount; iEdge++) {
                html << "<td class=centered>" << integerMatrix(iOrientedRead, iEdge);
            }
        }
        html << "</table>";
    }



    // Use a SVD to cluster the OrientedReadIds.

    // Create a matrix containing doubles.
    using Matrix = boost::numeric::ublas::matrix<double, boost::numeric::ublas::column_major>;
    Matrix matrix(orientedReadCount, edgeCount);
    for(uint64_t iOrientedRead=0; iOrientedRead<orientedReadCount; iOrientedRead++) {
        for(uint64_t iEdge=0; iEdge<edgeCount; iEdge++) {
            matrix(iOrientedRead, iEdge) = double(integerMatrix(iOrientedRead, iEdge));
        }
    }

    // Shift all the columns so they have zero average.
    for(uint64_t iEdge=0; iEdge<edgeCount; iEdge++) {
        double sum = 0.;
        for(uint64_t iOrientedRead=0; iOrientedRead<orientedReadCount; iOrientedRead++) {
            sum += matrix(iOrientedRead, iEdge);
        }
        const double average = sum / double(orientedReadCount);
        for(uint64_t iOrientedRead=0; iOrientedRead<orientedReadCount; iOrientedRead++) {
            matrix(iOrientedRead, iEdge) -= average;
        }
    }

    // Shift all the rows so they have zero average.
    for(uint64_t iOrientedRead=0; iOrientedRead<orientedReadCount; iOrientedRead++) {
        double sum = 0.;
        for(uint64_t iEdge=0; iEdge<edgeCount; iEdge++) {
            sum += matrix(iOrientedRead, iEdge);
        }
        const double average = sum / double(edgeCount);
        for(uint64_t iEdge=0; iEdge<edgeCount; iEdge++) {
            matrix(iOrientedRead, iEdge) -= average;
        }
    }

    // Compute the SVD.
    vector<double> singularValues;
    Matrix leftSingularVectors;
    Matrix rightSingularVectors;
    dgesvd(matrix, singularValues, leftSingularVectors, rightSingularVectors);



    // Write SVD results.
    if(html) {
        // The actual number of singular values/vectors we will write out.
        const uint64_t singularValueCount = min(6UL, singularValues.size());

        html << std::fixed << std::setprecision(3);

        // Singular values.
        html <<
            "<h3>Singular values</h3><table>";
        for(uint64_t j=0; j<singularValueCount; j++) {
            html << "<tr><th>S<sub>" << j << "</sub><td class=centered>" << singularValues[j];
        }
        html << "</table>";

        // Left singular vectors.
        html <<
            "<h3>Left singular vectors</h3><table>"
            "<tr><th>Oriented<br>read<br>id";
        for(uint64_t j=0; j<singularValueCount; j++) {
            html << "<th>L<sub>" << j << "</sub>";
        }
        for(uint64_t j=0; j<singularValueCount; j++) {
            html << "<th>S<sub>" << j << "</sub>L<sub>" << j << "</sub>";
        }
        html << "\n";
        for(uint64_t iOrientedRead=0; iOrientedRead<orientedReadCount; iOrientedRead++) {
            const OrientedReadId orientedReadId = orientedReadIds[iOrientedRead];
            html << "<tr><th>" << orientedReadId;

            // Without scaling.
            for(uint64_t j=0; j<singularValueCount; j++) {
                html << "<td class=centered>" << leftSingularVectors(iOrientedRead, j);
            }

            // With scaling.
            for(uint64_t j=0; j<singularValueCount; j++) {
                html << "<td class=centered>" << singularValues[j] * leftSingularVectors(iOrientedRead, j);
            }

            html << "\n";
        }
        html << "</table>";
    }
}
