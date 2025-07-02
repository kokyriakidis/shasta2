// Shasta.
#include "SearchGraph.hpp"
#include "transitiveReduction.hpp"
using namespace shasta;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>

// Standard library.
#include <fstream.hpp>



SearchGraph::SearchGraph(
    const AssemblyGraph& assemblyGraph,
    uint64_t minCoverage) :
    assemblyGraph(assemblyGraph)
{
    createVertices();
    createEdges(minCoverage);
    transitiveReductionAny(*this);

    cout << "The search graph has " << num_vertices(*this) <<
        " vertices and " << num_edges(*this) << " edges" << endl;
}


void SearchGraph::createVertices()
{
    SearchGraph& searchGraph = *this;

    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        const vertex_descriptor v = add_vertex(SearchGraphVertex(e), searchGraph);
        vertexMap.insert({e, v});
    }
}



void SearchGraph::createEdges(uint64_t minCoverage)
{
    SearchGraph& searchGraph = *this;

    // Use local searches to find forward edge pairs and backward edge pairs.
    using EdgePair = pair<edge_descriptor, edge_descriptor>;
    vector<EdgePair> forwardPairs;
    vector<EdgePair> backwardPairs;

    vector<edge_descriptor> localEdges;
    BGL_FORALL_EDGES(e0, assemblyGraph, AssemblyGraph) {
        assemblyGraph.forwardLocalSearch(e0, minCoverage, localEdges);
        for(const edge_descriptor e1: localEdges) {
            forwardPairs.push_back({e0, e1});
        }
        assemblyGraph.backwardLocalSearch(e0, minCoverage, localEdges);
        for(const edge_descriptor e1: localEdges) {
            backwardPairs.push_back({e1, e0});
        }
    }



    // Find edge pairs that appear in both directions.
    const AssemblyGraph::OrderById orderById(assemblyGraph);
    sort(forwardPairs.begin(), forwardPairs.end(), orderById);
    sort(backwardPairs.begin(), backwardPairs.end(), orderById);

    vector<EdgePair> bidirectionalPairs;
    std::set_intersection(
        forwardPairs.begin(), forwardPairs.end(),
        backwardPairs.begin(), backwardPairs.end(),
        back_inserter(bidirectionalPairs),
        orderById);



    // Each bidirectional pair generates an edge.
    for(const EdgePair& edgePair: bidirectionalPairs) {
        const AssemblyGraph::edge_descriptor e0 = edgePair.first;
        const AssemblyGraph::edge_descriptor e1 = edgePair.second;
        add_edge(vertexMap[e0], vertexMap[e1], searchGraph);
    }
}



void SearchGraph::writeGraphviz(const string& fileName) const
{
    ofstream s(fileName);
    writeGraphviz(s);
}



void SearchGraph::writeGraphviz(ostream& s) const
{
    const SearchGraph& searchGraph = *this;

    s << "digraph SearchGraph {\n";

    BGL_FORALL_VERTICES(v, searchGraph, SearchGraph) {
        const AssemblyGraph::edge_descriptor e = searchGraph[v].e;
        s << assemblyGraph[e].id << ";\n";
    }

    BGL_FORALL_EDGES(se, searchGraph, SearchGraph) {
        const vertex_descriptor sv0 = source(se, searchGraph);
        const vertex_descriptor sv1 = target(se, searchGraph);
        const AssemblyGraph::edge_descriptor e0 = searchGraph[sv0].e;
        const AssemblyGraph::edge_descriptor e1 = searchGraph[sv1].e;
        s << assemblyGraph[e0].id << "->" << assemblyGraph[e1].id << ";\n";
    }

    s << "}\n";
}
