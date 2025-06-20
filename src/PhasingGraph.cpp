// Shasta.
#include "PhasingGraph.hpp"
#include "color.hpp"
#include "orderVectors.hpp"
#include "SHASTA_ASSERT.hpp"
using namespace shasta;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>
#include <boost/pending/disjoint_sets.hpp>

// Standard library.
#include <algorithm.hpp>
#include <fstream.hpp>



void PhasingGraph::addVertex(uint64_t position)
{
    // Make sure the vertexTable can store this position.
    if(position >= vertexTable.size()) {
        vertexTable.resize(position+1, null_vertex());
    }

    // Check that we don't already have a vertex at this position.
    SHASTA_ASSERT(vertexTable[position] == null_vertex());

    // Add the vertex.
    const vertex_descriptor v = boost::add_vertex(PhasingGraphVertex(position), *this);
    vertexTable[position] = v;
}



void PhasingGraph::addEdge(
    uint64_t position0,
    uint64_t position1,
    const TangleMatrix::Hypothesis& bestHypothesis)
{
    SHASTA_ASSERT(position0 < vertexTable.size());
    const vertex_descriptor v0 = vertexTable[position0];
    SHASTA_ASSERT(v0 != null_vertex());

    SHASTA_ASSERT(position1 < vertexTable.size());
    const vertex_descriptor v1 = vertexTable[position1];
    SHASTA_ASSERT(v1 != null_vertex());

    boost::add_edge(v0, v1, PhasingGraphEdge(bestHypothesis), *this);
}



void PhasingGraph::writeGraphviz(const string& fileName) const
{
    const PhasingGraph& phasingGraph = *this;

    ofstream dot(fileName);
    dot << "digraph PhasingGraph {\n";

    BGL_FORALL_VERTICES(v, phasingGraph, PhasingGraph) {
        const PhasingGraphVertex& vertex = phasingGraph[v];
        const string color = hslToRgbString(double(vertex.componentId) / double(components.size()), 0.75, 0.6);
        dot << vertex.position << " [label=\"" <<
            vertex.position << "\\n" <<
            vertex.componentId << "\\n" <<
            vertex.pathLength << "\""
            " style=filled fillcolor=\"" << color << "\"];\n";
    }

    BGL_FORALL_EDGES(e, phasingGraph, PhasingGraph) {
        const vertex_descriptor v0 = source(e, phasingGraph);
        const vertex_descriptor v1 = target(e, phasingGraph);

        dot << phasingGraph[v0].position << "->";
        dot << phasingGraph[v1].position;
        if(phasingGraph[e].isShortestPathEdge) {
            dot << " [color=DarkOrange]";
        }
        dot << ";\n";
    }

    dot << "}\n";
}



// Remove isolated vertices and return the number of such vertices that were removed.
uint64_t PhasingGraph::removeIsolatedVertices()
{
    PhasingGraph& phasingGraph = *this;

    vector<vertex_descriptor> verticesToBeRemoved;
    BGL_FORALL_VERTICES(v, phasingGraph, PhasingGraph) {
        if((in_degree(v, phasingGraph) == 0) and (out_degree(v, phasingGraph)== 0)) {
            verticesToBeRemoved.push_back(v);
        }
    }

    for(const vertex_descriptor v: verticesToBeRemoved) {
        vertexTable[phasingGraph[v].position] = null_vertex();
        boost::remove_vertex(v, phasingGraph);
    }

    return verticesToBeRemoved.size();
}



// Compute connected components consisting of at least two vertices.
// Each connected component is a sorted vector of positions in the SuperbubbleChain.
// They are returned sorted by decreasing size.
void PhasingGraph::computeConnectedComponents()
{
    PhasingGraph& phasingGraph = *this;
    const uint64_t n = vertexTable.size();

    // Initialize the disjoint sets data structure.
    vector<uint64_t> rank(n);
    vector<uint64_t> parent(n);
    boost::disjoint_sets<uint64_t*, uint64_t*> disjointSets(&rank[0], &parent[0]);
    for(uint64_t position=0; position<n; position++) {
        disjointSets.make_set(position);
    }

    // Compute the disjoint sets.
    BGL_FORALL_EDGES(e, phasingGraph, PhasingGraph) {
        const vertex_descriptor v0 = source(e, phasingGraph);
        const vertex_descriptor v1 = target(e, phasingGraph);

        const uint64_t position0 = phasingGraph[v0].position;
        const uint64_t position1 = phasingGraph[v1].position;

        SHASTA_ASSERT(position0 < n);
        SHASTA_ASSERT(position1 < n);

        disjointSets.union_set(position0, position1);
    }

    // Gather the vertices in each connected component.
    vector< vector<uint64_t> > componentTable(n);
    for(uint64_t position=0; position<n; position++) {
        const uint64_t componentId = disjointSets.find_set(position);
        componentTable[componentId].push_back(position);
    }

    components.clear();
    for(const vector<uint64_t>& component: componentTable) {
        if(component.size() > 1) {
            components.push_back(component);
        }
    }

    // Sort them by decreasing size.
    sort(components.begin(), components.end(), OrderVectorsByDecreasingSize<uint64_t>());

    // Store in each vertex the component it belongs to.
    for(uint64_t componentId=0; componentId<components.size(); componentId++) {
        const vector<uint64_t> component = components[componentId];
        for(const uint64_t position: component) {
            const vertex_descriptor v = vertexTable[position];
            SHASTA_ASSERT(v != null_vertex());
            phasingGraph[v].componentId = componentId;
        }
    }

}



// Find the longest path in each connected component.
// We exploit the fact that the PhasingGraph is acyclic and
// topologically sorted by constuction.
// https://en.wikipedia.org/wiki/Longest_path_problem#Acyclic_graphs
void PhasingGraph::findLongestPaths()
{
    PhasingGraph& phasingGraph = *this;

    // For all vertices, compute the maximum path length from a source vertex.
    for(uint64_t position0=0; position0<vertexTable.size(); position0++) {
        const vertex_descriptor v0 = vertexTable[position0];
        if(v0 == null_vertex()) {
            continue;
        }

        PhasingGraphVertex& vertex0 = phasingGraph[v0];
        if(in_degree(v0, phasingGraph) == 0) {
            vertex0.pathLength = 0;
        }

        const uint64_t pathLength0 = vertex0.pathLength;
        const uint64_t pathLength1 = pathLength0 + 1;

        BGL_FORALL_OUTEDGES(v0, e, phasingGraph, PhasingGraph) {
            const vertex_descriptor v1 = target(e, phasingGraph);
            PhasingGraphVertex& vertex1 = phasingGraph[v1];

            if(vertex1.pathLength == invalid<uint64_t>) {
                vertex1.pathLength = pathLength1;
            } else {
                vertex1.pathLength = max(vertex1.pathLength, pathLength1);
            }

        }
    }



    // Now we can compute the longest path in each connected component.
    longestPaths.resize(components.size());
    for(uint64_t componentId=0; componentId<components.size(); componentId++) {
        const vector<uint64_t>& component = components[componentId];
        vector<edge_descriptor>& path = longestPaths[componentId];

        // Find a vertex with the largest pathLength. This will be
        // the last vertex of our path.
        uint64_t maxPathLength = 0;
        vertex_descriptor vLast = null_vertex();
        for(const uint64_t position: component) {
            const vertex_descriptor v = vertexTable[position];
            const PhasingGraphVertex& vertex = phasingGraph[v];
            SHASTA_ASSERT(v != null_vertex());

            if(vLast == null_vertex()) {
                vLast = v;
                maxPathLength = vertex.pathLength;
            } else {
                if(vertex.pathLength > maxPathLength) {
                    vLast = v;
                    maxPathLength = vertex.pathLength;
                }
            }
        }

        // Now we can construct the path walking back from vLast.
        vertex_descriptor v1 = vLast;
        while(in_degree(v1, phasingGraph) > 0) {
            const PhasingGraphVertex& vertex1 = phasingGraph[v1];
            const uint64_t pathLength1 = vertex1.pathLength;
            const uint64_t pathLength0 = pathLength1 - 1;
            bool found = false;
            BGL_FORALL_INEDGES(v1, e, phasingGraph, PhasingGraph) {
                const vertex_descriptor v0 = source(e, phasingGraph);
                if(phasingGraph[v0].pathLength == pathLength0) {
                    phasingGraph[e].isShortestPathEdge = true;
                    path.push_back(e);
                    v1 = v0;
                    found = true;
                    break;
                }
            }
            SHASTA_ASSERT(found);
        }
        std::reverse(path.begin(), path.end());

    }

}

