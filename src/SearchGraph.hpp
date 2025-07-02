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
    SearchGraphVertex() {}
};



class shasta::SearchGraph : public SearchGraphBaseClass {
public:

    // Initial construction from the AssemblyGraph.
    SearchGraph(
        const AssemblyGraph&,
        uint64_t minCoverage);

    // Construction from connected component of the SearchGraph.
    SearchGraph(const SearchGraph&, const vector<vertex_descriptor>&);

    void writeGraphviz(const string& fileName) const;
    void writeGraphviz(ostream&) const;

    // This computes connected components and creates a new SearchGraph
    // for each non-trivial connected component.
    void computeConnectedComponents(vector<SearchGraph>&) const;

    // This removes out-edges of vertices with out_degree > 1
    // and in-edges of vertices with in-degree > 1.
    // After this operation, the SearchGraph consists of a set
    // of linear chains.
    void removeBranches();

private:
    const AssemblyGraph& assemblyGraph;
    std::map<AssemblyGraph::edge_descriptor, vertex_descriptor> vertexMap;

    void createVertices();
    void createEdges(uint64_t minCoverage);

};
