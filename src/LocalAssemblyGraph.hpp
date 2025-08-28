#pragma once

// Shasta.
#include "AssemblyGraph.hpp"

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>

namespace shasta {
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


class shasta::LocalAssemblyGraphVertex {
public:
    AssemblyGraph::vertex_descriptor v;
    uint64_t distance = invalid<uint64_t>;

    // Fields used by approximateTopologicalSort.
    uint64_t color = invalid<uint64_t>;
    uint64_t rank = invalid<uint64_t>;
};



class shasta::LocalAssemblyGraphEdge {
public:
    AssemblyGraph::edge_descriptor e;

    // The default constructor is only needed to avoid a warning.
    LocalAssemblyGraphEdge() : e(0, 0, 0) {}

    // Field used by approximateTopologicalSort.
    bool isDagEdge = false;
};



class shasta::LocalAssemblyGraph : public LocalAssemblyGraphBaseClass {
public:
    LocalAssemblyGraph(
        const AssemblyGraph&,
        vector<AssemblyGraph::vertex_descriptor> startVertices,
        uint64_t maxDistance);
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
};
