// Shasta2.
#include "TangleGraph.hpp"
#include "AssemblyGraph.hpp"
#include "html.hpp"
#include "TangleMatrix.hpp"
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



    // Fill in the reverse complement of each vertex.
    BGL_FORALL_VERTICES(tv, tangleGraph, TangleGraph) {
        const AssemblyGraph::vertex_descriptor av = tangleGraph[tv].assemblyGraphVertices.front();
        const AssemblyGraph::vertex_descriptor avRc = assemblyGraph[av].vRc;
        tangleGraph[tv].vRc = vertexMap.at(avRc);
    }

    // Sanity check of the reverse complement of each vertex.
    BGL_FORALL_VERTICES(v, tangleGraph, TangleGraph) {
        const TangleGraphVertex& vertex = tangleGraph[v];
        const vertex_descriptor vRc = vertex.vRc;
        const TangleGraphVertex& vertexRc = tangleGraph[vRc];
        SHASTA2_ASSERT(vertexRc.vRc == v);
    }


    writeVertices("TangleGraphVertices.csv");
    writeEdges("TangleGraphEdges.csv");
    writeTangleMatrices("TangleMatrices.html");
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



void TangleGraph::writeVertices(const string& fileName) const
{
    ofstream csv(fileName);
    writeVertices(csv);
}



void TangleGraph::writeVertices(ostream& csv) const
{
    const TangleGraph& tangleGraph = *this;

    csv << "Id,IdRc,VertexCount,Entrances,Exits,\n";

    BGL_FORALL_VERTICES(tv, tangleGraph, TangleGraph) {
        const TangleGraphVertex& vertex = tangleGraph[tv];
        vector<AssemblyGraph::edge_descriptor> entrances;
        vector<AssemblyGraph::edge_descriptor> exits;
        getEntrances(tv, entrances);
        getExits(tv, exits);

        csv <<
            vertex.id << "," <<
            tangleGraph[vertex.vRc].id << "," <<
            vertex.assemblyGraphVertices.size() << ",";

        for(uint64_t i=0; i<entrances.size(); i++) {
            const AssemblyGraph::edge_descriptor entrance = entrances[i];
            if(i != 0) {
                csv << " ";
            }
            csv << assemblyGraph[entrance].id;
        }
        csv << ",";

        for(uint64_t i=0; i<exits.size(); i++) {
            const AssemblyGraph::edge_descriptor exit = exits[i];
            if(i != 0) {
                csv << " ";
            }
            csv << assemblyGraph[exit].id;
        }
        csv << ",";

        for(const AssemblyGraph::vertex_descriptor av: vertex.assemblyGraphVertices) {
            csv << assemblyGraph[av].id << ",";
        }
        csv << "\n";
    }

}



void TangleGraph::writeEdges(const string& fileName) const
{
    ofstream csv(fileName);
    writeEdges(csv);
}



void TangleGraph::writeEdges(ostream& csv) const
{
    const TangleGraph& tangleGraph = *this;

    csv << "Id,Length,\n";

    BGL_FORALL_EDGES(ev, tangleGraph, TangleGraph) {
        const TangleGraphEdge& edge = tangleGraph[ev];
        const BubbleChain& bubbleChain = edge.bubbleChain;
        csv <<
            edge.id << "," <<
            bubbleChain.maxLength(assemblyGraph) << ",";
        for(const Bubble& bubble: bubbleChain) {
            for(uint64_t i=0; i<bubble.edges.size(); i++) {
                const AssemblyGraph::edge_descriptor e = bubble.edges[i];
                if(i != 0) {
                    csv << " ";
                }
                csv << assemblyGraph[e].id;
            }
            csv << ",";
        }
        csv << "\n";
    }
}


void TangleGraph::getEntrances(
    vertex_descriptor v,
    vector<AssemblyGraph::edge_descriptor>& entrances) const
{
    const TangleGraph& tangleGraph = *this;

    entrances.clear();
    BGL_FORALL_INEDGES(v, e, tangleGraph, TangleGraph) {
        const BubbleChain& bubbleChain = tangleGraph[e].bubbleChain;
        const Bubble& bubble = bubbleChain.back();
        std::ranges::copy(bubble.edges, back_inserter(entrances));
    }

    sort(entrances.begin(), entrances.end(), assemblyGraph.orderById);
}



void TangleGraph::getExits(
    vertex_descriptor v,
    vector<AssemblyGraph::edge_descriptor>& exits) const
{
    const TangleGraph& tangleGraph = *this;

    exits.clear();
    BGL_FORALL_OUTEDGES(v, e, tangleGraph, TangleGraph) {
        const BubbleChain& bubbleChain = tangleGraph[e].bubbleChain;
        const Bubble& bubble = bubbleChain.front();
        std::ranges::copy(bubble.edges, back_inserter(exits));
    }

    sort(exits.begin(), exits.end(), assemblyGraph.orderById);
}



void TangleGraph::writeTangleMatrices(const string& fileName) const
{
    ofstream html(fileName);
    writeTangleMatrices(html);
}



void TangleGraph::writeTangleMatrices(ostream& html) const
{
    const TangleGraph& tangleGraph = *this;
    ostream noOutput(0);
    writeHtmlBegin(html, "Tangle matrices");

    BGL_FORALL_VERTICES(v, tangleGraph, TangleGraph) {
        vector<AssemblyGraph::edge_descriptor> entrances;
        vector<AssemblyGraph::edge_descriptor> exits;
        getEntrances(v, entrances);
        getExits(v, exits);
        const TangleMatrix tangleMatrix(assemblyGraph, entrances, exits, noOutput);

        html << "<h2>Tangle graph vertex " << tangleGraph[v].id << "</h2>" << endl;

        html << "<table><tr><th>";
        for(uint64_t j=0; j<exits.size(); j++) {
            html << "<th>" << assemblyGraph[exits[j]].id;
        }

        for(uint64_t i=0; i<entrances.size(); i++) {
            html << "<tr><th>" << assemblyGraph[entrances[i]].id;
            for(uint64_t j=0; j<exits.size(); j++) {
                html << "<td class=centered>" << tangleMatrix.tangleMatrix[i][j];
            }
        }

        html << "</table>";
    }

    writeHtmlEnd(html);

}
