#ifndef SHASTA_DOMINATOR_TREE_HPP
#define SHASTA_DOMINATOR_TREE_HPP

/*******************************************************************************


This provides a bug fix for the 3-argument version of
boost::lengauer_tarjan_dominator_tree.
The fixed version is in the shasta namespace.

The bug fix is as follows.

The original version contains

    std::vector<VerticesSizeType> dfnum(numOfVertices, 0);

The fixed version below has instead

    std::vector<VerticesSizeType> dfnum(numOfVertices, std::numeric_limits<VerticesSizeType>::max());

As a result of the bug, the original version produces incorrect results
if the graph contains vertices that are unreachable from the entrance.

lengauer_tarjan_dominator_tree_without_dfs, which is eventually called,
includes the following comment about initializing dfnum/dfnumMap:

   * @pre Unreachable nodes must be masked as
   *      (std::numeric_limits<VerticesSizeType>::max)() in dfnumMap.

The above fix implements this requirement.
The fixed version below works correctly even if the graph contains unreachable vertices.

*******************************************************************************/

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/graph/dominator_tree.hpp>

namespace shasta {

    template<class Graph, class DomTreePredMap>
    void lengauer_tarjan_dominator_tree(const Graph &g,
        const typename boost::graph_traits<Graph>::vertex_descriptor &entry,
        DomTreePredMap domTreePredMap)
    {
        using namespace boost;

        // typedefs
        typedef typename graph_traits<Graph>::vertex_descriptor Vertex;
        typedef typename graph_traits<Graph>::vertices_size_type VerticesSizeType;
        typedef typename property_map<Graph, vertex_index_t>::const_type IndexMap;
        typedef iterator_property_map<
            typename std::vector<VerticesSizeType>::iterator, IndexMap> TimeMap;
        typedef iterator_property_map<typename std::vector<Vertex>::iterator,
            IndexMap> PredMap;

        // Make property maps
        const VerticesSizeType numOfVertices = num_vertices(g);
        if (numOfVertices == 0)
            return;

        const IndexMap indexMap = get(vertex_index, g);

        // ********* THE BUG FIX IS HERE **********************
        // The original versionwas filling the vector with 0.
        std::vector<VerticesSizeType> dfnum(numOfVertices,
            std::numeric_limits < VerticesSizeType > ::max());

        TimeMap dfnumMap(make_iterator_property_map(dfnum.begin(), indexMap));

        std::vector<Vertex> parent(numOfVertices,
            graph_traits < Graph > ::null_vertex());
        PredMap parentMap(make_iterator_property_map(parent.begin(), indexMap));

        std::vector<Vertex> verticesByDFNum(parent);

        // Run main algorithm
        boost::lengauer_tarjan_dominator_tree(g, entry, indexMap, dfnumMap,
            parentMap, verticesByDFNum, domTreePredMap);
    }



    // The above requires a Graph that uses vecS to store the vertices.
    // For graphs that don't do that, I was not able to get the 7-argument
    // version of boost::lengauer_tarjan_dominator_tree to work.
    // As a workaround, the code below creates an AuxiliaryGraph
    // that uses vecS, computes the dominator tree for the AuxiliaryGraph,
    // then stores the result in the dominator field of the Graph vertex (required).
    template<class Graph>
    void lengauer_tarjan_dominator_tree_general(
        Graph &graph,
        typename boost::graph_traits<Graph>::vertex_descriptor entry)
    {
        using namespace boost;

        // Vertex descriptor for the original Graph.
        using V = typename Graph::vertex_descriptor;

        // Base class for AuxiliaryGraph.
        class AuxiliaryGraphVertex;
        using AuxiliaryGraph = boost::adjacency_list<
            boost::listS,
            boost::vecS,
            boost::bidirectionalS,
            AuxiliaryGraphVertex>;

        // Vertex descriptor for the AuxiliaryGraph.
        using AV = typename AuxiliaryGraph::vertex_descriptor;

        // Each vertex stores the corresponding vertex descriptor of the Graph.
        class AuxiliaryGraphVertex {
        public:
            V v;
            AV dominator = AuxiliaryGraph::null_vertex();
            AuxiliaryGraphVertex(V v = Graph::null_vertex()) : v(v) {}
         };

        // Create the AuxiliaryGraph.
        AuxiliaryGraph auxiliaryGraph;

        // Create an AuxiliaryGraph vertex for each vertex of the Graph.
        std::map<V, AV> vertexMap;
        BGL_FORALL_VERTICES_T(v, graph, Graph) {
            const AV av = boost::add_vertex(AuxiliaryGraphVertex(v), auxiliaryGraph);
            vertexMap.insert(make_pair(v, av));
        }

        // Create an AuxiliaryGraph edge for each edge of the Graph.
        BGL_FORALL_EDGES_T(e, graph, Graph) {
            const V v0 = source(e, graph);
            const V v1 = target(e, graph);
            const AV av0 = vertexMap[v0];
            const AV av1 = vertexMap[v1];
            boost::add_edge(av0, av1, auxiliaryGraph);
        }

        // Compute the dominator tree of the AuxiliaryGraph.
        const auto itEntry = vertexMap.find(entry);
        SHASTA_ASSERT(itEntry != vertexMap.end());
        const AV auxEntry = itEntry->second;
        shasta::lengauer_tarjan_dominator_tree(
            auxiliaryGraph,
            auxEntry,
            boost::get(&AuxiliaryGraphVertex::dominator, auxiliaryGraph));

        // Store it in the dominator field of the vertices of the original graph.
        BGL_FORALL_VERTICES_T(v, graph, Graph) {
            const AV av = vertexMap[v];
            const AV auxDominator = auxiliaryGraph[av].dominator;
            if(auxDominator == AuxiliaryGraph::null_vertex()){
                graph[v].dominator = Graph::null_vertex();
            } else {
                graph[v].dominator = auxiliaryGraph[auxDominator].v;
            }
        }
    }
}

#endif

