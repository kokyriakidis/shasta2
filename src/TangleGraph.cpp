// Shasta2.
#include "TangleGraph.hpp"
#include "AssemblyGraph.hpp"
using namespace shasta2;

// Standard library.
#include "fstream.hpp"
#include "iostream.hpp"



TangleGraph::TangleGraph(
    const AssemblyGraph& assemblyGraph,
    const vector<BubbleChain>& bubbleChains,
    const vector< vector<vertex_descriptor> >& tangles) :
    assemblyGraph(assemblyGraph)
{
    TangleGraph& tangleGraph = *this;

    // Generate a vertex for each tangle and keep track of the
    // AssemblyGraph vertices in each tangle.
    uint64_t nextVertexId = 0;
    std::map<AssemblyGraph::vertex_descriptor, TangleGraph::vertex_descriptor> vertexMap;
    for(const vector<vertex_descriptor>& tangle: tangles) {
        const vertex_descriptor tv = boost::add_vertex(TangleGraphVertex(nextVertexId++, tangle), tangleGraph);
        for(const AssemblyGraph::vertex_descriptor v: tangle) {
            vertexMap.insert(make_pair(v, tv));
        }
    }


    // Now generate an edge for each BubbleChain.
    for(uint64_t bubbleChainId=0; bubbleChainId<bubbleChains.size(); bubbleChainId++) {
        const BubbleChain& bubbleChain = bubbleChains[bubbleChainId];
        const AssemblyGraph::vertex_descriptor v0 = bubbleChain.front().v0;
        const AssemblyGraph::vertex_descriptor v1 = bubbleChain.back().v1;

        TangleGraph::vertex_descriptor tv0;
        auto it0 = vertexMap.find(v0);
        if(it0 == vertexMap.end()) {
            tv0 = boost::add_vertex(TangleGraphVertex(nextVertexId++, v0), tangleGraph);
            vertexMap.insert(make_pair(v0, tv0));
        } else {
            tv0 = it0->second;
        }

        TangleGraph::vertex_descriptor tv1;
        auto it1 = vertexMap.find(v1);
        if(it1 == vertexMap.end()) {
            tv1 = boost::add_vertex(TangleGraphVertex(nextVertexId++, v1), tangleGraph);
            vertexMap.insert(make_pair(v1, tv1));
        } else {
            tv1 = it1->second;
        }

       boost::add_edge(tv0, tv1, TangleGraphEdge(bubbleChainId, bubbleChain), tangleGraph);
    }

    writeGraphviz("TangleGraph.dot");
    cout << "The TangleGraph has " << num_vertices(tangleGraph) <<
        " vertices and " << num_edges(tangleGraph) << " edges." << endl;
}



void TangleGraph::writeGraphviz(const string& fileName) const
{
    ofstream dot(fileName);
    writeGraphviz(dot);
}



void TangleGraph::writeGraphviz(ostream& dot) const
{
    const TangleGraph& tangleGraph = *this;

    dot << "digraph TangleGraph {\n";

    BGL_FORALL_VERTICES(v, tangleGraph, TangleGraph) {
        dot << tangleGraph[v].id << ";\n";
    }

    BGL_FORALL_EDGES(e, tangleGraph, TangleGraph) {
        const vertex_descriptor v0 = source(e, tangleGraph);
        const vertex_descriptor v1 = target(e, tangleGraph);
        dot <<
            tangleGraph[v0].id << "->" <<
            tangleGraph[v1].id <<
            " [label=\"" << tangleGraph[e].id << "\\n" <<
            tangleGraph[e].bubbleChain.maxLength(assemblyGraph) << "\"]"
            ";\n";
    }

    dot << "}\n";
}

