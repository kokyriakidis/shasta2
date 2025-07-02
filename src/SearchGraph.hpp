#pragma once

#include "AssemblyGraph.hpp"

namespace shasta {
    class SearchGraphVertex;
    class SearchGraph;

    using SearchGraphBaseClass = boost::adjacency_list<
        boost::listS,
        boost::listS,
        boost::bidirectionalS,
        SearchGraphVertex>;
}



class shasta::SearchGraphVertex {
public:
    AssemblyGraph::edge_descriptor e;

    SearchGraphVertex(AssemblyGraph::edge_descriptor e) : e(e) {}
};



class shasta::SearchGraph : public SearchGraphBaseClass {
public:

    SearchGraph(
        const AssemblyGraph&,
        uint64_t minCoverage);

    void writeGraphviz(const string& fileName) const;
    void writeGraphviz(ostream&) const;

private:
    const AssemblyGraph& assemblyGraph;
    std::map<AssemblyGraph::edge_descriptor, vertex_descriptor> vertexMap;

    void createVertices();
    void createEdges(uint64_t minCoverage);

};
