// Shasta.
#include "SearchGraph.hpp"
#include "deduplicate.hpp"
using namespace shasta;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>
#include <boost/pending/disjoint_sets.hpp>

// Standard library.
#include <fstream.hpp>



SearchGraph::SearchGraph(
    const AssemblyGraph& assemblyGraph,
    uint64_t lowCoverageThreshold,
    uint64_t highCoverageThreshold) :
    assemblyGraph(assemblyGraph)
{
    createVertices();
    createEdges(lowCoverageThreshold, highCoverageThreshold);

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



void SearchGraph::createEdges(
    uint64_t lowCoverageThreshold,
    uint64_t highCoverageThreshold)
{
    SearchGraph& searchGraph = *this;

    // Use local searches to find forward edge pairs and backward edge pairs.
    using EdgePair = pair<edge_descriptor, edge_descriptor>;
    vector<EdgePair> forwardPairs;
    vector<EdgePair> backwardPairs;

    vector<edge_descriptor> localEdges;
    BGL_FORALL_EDGES(e0, assemblyGraph, AssemblyGraph) {
        assemblyGraph.forwardLocalSearch(e0, lowCoverageThreshold, highCoverageThreshold, localEdges);
        for(const edge_descriptor e1: localEdges) {
            forwardPairs.push_back({e0, e1});
        }
        assemblyGraph.backwardLocalSearch(e0, lowCoverageThreshold, highCoverageThreshold, localEdges);
        for(const edge_descriptor e1: localEdges) {
            backwardPairs.push_back({e1, e0});
        }
    }



    // Find edge pairs that appear in both directions.
    sort(forwardPairs.begin(), forwardPairs.end(), assemblyGraph.orderById);
    sort(backwardPairs.begin(), backwardPairs.end(), assemblyGraph.orderById);

    vector<EdgePair> bidirectionalPairs;
    std::set_intersection(
        forwardPairs.begin(), forwardPairs.end(),
        backwardPairs.begin(), backwardPairs.end(),
        back_inserter(bidirectionalPairs),
        assemblyGraph.orderById);



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



// This computes connected components and creates a new SearchGraph
// for each non-trivial connected component.
void SearchGraph::computeConnectedComponents(vector<SearchGraph>& components) const
{
    const SearchGraph& searchGraph  = *this;

    // Map vertices to integer.
    std::map<vertex_descriptor, uint64_t> indexMap;
    uint64_t vertexIndex = 0;
    BGL_FORALL_VERTICES(v, searchGraph, SearchGraph) {
        indexMap.insert({v, vertexIndex++});
    }
    const uint64_t vertexCount = indexMap.size();

    // Initialize the disjoint sets.
    vector<uint64_t> rank(vertexCount);
    vector<uint64_t> parent(vertexCount);
    boost::disjoint_sets<uint64_t*, uint64_t*> disjointSets(&rank[0], &parent[0]);
    for(uint64_t i=0; i<vertexCount; i++) {
        disjointSets.make_set(i);
    }

    // Compute connected components.
    BGL_FORALL_EDGES(e, searchGraph, SearchGraph) {
        const vertex_descriptor v0 = source(e, searchGraph);
        const vertex_descriptor v1 = target(e, searchGraph);
        const uint64_t index0 = indexMap[v0];
        const uint64_t index1 = indexMap[v1];
        disjointSets.union_set(index0, index1);
    }

    // Gather the vertices in each connected component.
    vector < vector<vertex_descriptor> > componentTable(vertexCount);
    BGL_FORALL_VERTICES(v, searchGraph, SearchGraph) {
        const uint64_t index = indexMap[v];
        const uint64_t componentId = disjointSets.find_set(index);
        componentTable[componentId].push_back(v);
    }



    // Create and store a new SearchGraph for each non-trivial connected components.
    components.clear();
    for(const vector<vertex_descriptor>& componentVertices: componentTable) {
        if(componentVertices.size() > 1) {
            components.emplace_back(*this, componentVertices);
        }
    }
}


// Construction from connected component of the SearchGraph.
SearchGraph::SearchGraph(
    const SearchGraph& searchGraph,
    const vector<vertex_descriptor>& componentVertices) :
    assemblyGraph(searchGraph.assemblyGraph)
{
    SearchGraph& component = *this;

    // Add the vertices.
    for(const vertex_descriptor vGlobal: componentVertices) {
        const AssemblyGraph::edge_descriptor e = searchGraph[vGlobal].e;
        const vertex_descriptor v = add_vertex(SearchGraphVertex(e), component);
        vertexMap.insert({e, v});
    }

    // Add the edges.
    for(const vertex_descriptor v0Global: componentVertices) {
        const AssemblyGraph::edge_descriptor e0 = searchGraph[v0Global].e;
        const auto it0 = vertexMap.find(e0);
        SHASTA2_ASSERT(it0 != vertexMap.end());
        const vertex_descriptor v0 = it0->second;
        BGL_FORALL_OUTEDGES(v0Global, eGlobal, searchGraph, SearchGraph) {
            const vertex_descriptor v1Global = target(eGlobal, searchGraph);
            const AssemblyGraph::edge_descriptor e1 = searchGraph[v1Global].e;
            const auto it1 = vertexMap.find(e1);
            SHASTA2_ASSERT(it1 != vertexMap.end());
            const vertex_descriptor v1 = it1->second;
            add_edge(v0, v1, component);
        }
    }


}




// This removes out-edges of vertices with out_degree > 1
// and in-edges of vertices with in-degree > 1.
// After this operation, the SearchGraph consists of a set
// of linear chains.
void SearchGraph::removeBranches()
{
    SearchGraph& searchGraph = *this;

    vector<edge_descriptor> edgesToBeRemoved;
    BGL_FORALL_VERTICES(v, searchGraph, SearchGraph) {
        if(out_degree(v, searchGraph) > 1) {
            BGL_FORALL_OUTEDGES(v, e, searchGraph, SearchGraph) {
                edgesToBeRemoved.push_back(e);
            }
        }
        if(in_degree(v, searchGraph) > 1) {
            BGL_FORALL_INEDGES(v, e, searchGraph, SearchGraph) {
                edgesToBeRemoved.push_back(e);
            }
        }
    }
    deduplicate(edgesToBeRemoved);

    for(const edge_descriptor e: edgesToBeRemoved) {
        boost::remove_edge(e, searchGraph);
    }
}
