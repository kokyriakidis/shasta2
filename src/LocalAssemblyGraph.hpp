#pragma once

// Shasta.
#include "AssemblyGraph.hpp"

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>

namespace shasta2 {
    class LocalAssemblyGraph;
    class LocalAssemblyGraphVertex;
    class LocalAssemblyGraphEdge;

    using LocalAssemblyGraphBaseClass = boost::adjacency_list<
        boost::listS,
        boost::listS,
        boost::bidirectionalS,
        LocalAssemblyGraphVertex,
        LocalAssemblyGraphEdge>;
}


class shasta2::LocalAssemblyGraphVertex {
public:
    AssemblyGraph::vertex_descriptor v;

    LocalAssemblyGraphVertex(AssemblyGraph::vertex_descriptor v) : v(v) {}
    LocalAssemblyGraphVertex() : v(LocalAssemblyGraphBaseClass::null_vertex()) {}

    uint64_t distance = invalid<uint64_t>;

    // Fields used by approximateTopologicalSort.
    uint64_t color = invalid<uint64_t>;
    uint64_t rank = invalid<uint64_t>;
};



class shasta2::LocalAssemblyGraphEdge {
public:
    AssemblyGraph::edge_descriptor e;

    LocalAssemblyGraphEdge(AssemblyGraph::edge_descriptor e) : e(e) {}
    LocalAssemblyGraphEdge() : e(0, 0, 0) {}

    // Field used by approximateTopologicalSort.
    bool isDagEdge = false;
};



class shasta2::LocalAssemblyGraph : public LocalAssemblyGraphBaseClass {
public:

    // Constructor from a set of AssemblyGraph start vertices and maximum distance.
    LocalAssemblyGraph(
        const AssemblyGraph&,
        vector<AssemblyGraph::vertex_descriptor> startVertices,
        uint64_t maxDistance);

    // Constructor from a set of AssemblyGraph edges.
    LocalAssemblyGraph(
        const AssemblyGraph&,
        vector<AssemblyGraph::edge_descriptor> edges);

    // Default constructor.
    LocalAssemblyGraph() {}

    void writeHtml(
        ostream& html,
        const AssemblyGraph&,
        uint64_t maxDistance) const;
    void writeGraphviz(
        const string& fileName,
        const AssemblyGraph&,
        uint64_t maxDistance) const;
    void writeGraphviz(
        ostream&,
        const AssemblyGraph&,
        uint64_t maxDistance) const;

    void approximateTopologicalSort(const AssemblyGraph&);
    vector<vertex_descriptor> verticesByRank;

    void analyze(ostream& html, const AssemblyGraph&) const;
};
