// Shasta.
#include "TransitionGraph.hpp"
#include "deduplicate.hpp"
#include "performanceLog.hpp"
#include "timestamp.hpp"
using namespace shasta;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>

// Standard library.
#include <fstream.hpp>



TransitionGraph::TransitionGraph(
    const Anchors& anchors,
    const AnchorGraph& anchorGraph)
{
    const uint64_t minTransitionGraphEdgeCoverage = 10;
    TransitionGraph& transitionGraph = *this;

    // Create the vertices, one for each AnchorGraph edge.
    uint64_t vertexId = 0;;
    BGL_FORALL_EDGES(eAnchorGraph, anchorGraph, AnchorGraph) {
        const vertex_descriptor v = add_vertex(TransitionGraphVertex(eAnchorGraph, vertexId++), transitionGraph);
        vertexMap.insert(make_pair(eAnchorGraph, v));
    }

    // Compute AnchorGraph edge journeys.
    // The edge journey of an OrientedReadId is the sequence of
    // AnchorGraph edges visited by the OrientedReadId.
    vector< vector<AnchorGraph::edge_descriptor> > edgeJourneys;
    performanceLog << timestamp << "AnchorGraph edge journey creation begins." << endl;
    anchorGraph.computeEdgeJourneys(anchors, edgeJourneys);
    performanceLog << timestamp << "AnchorGraph edge journey creation ends." << endl;



    // Gather pairs of consecutive AnchorGraph edges visited by all OrientedReadIds.
    performanceLog << timestamp << "Finding edge journey pairs begins." << endl;
    vector< pair<AnchorGraph::edge_descriptor, AnchorGraph::edge_descriptor> > edgeJourneyPairs;
    for(const vector<AnchorGraph::edge_descriptor>& edgeJourney: edgeJourneys) {
        for(uint64_t i1=1; i1<edgeJourney.size(); i1++) {
            const uint64_t i0 = i1 - 1;
            const AnchorGraph::edge_descriptor e0 = edgeJourney[i0];
            const AnchorGraph::edge_descriptor e1 = edgeJourney[i1];
            edgeJourneyPairs.push_back(make_pair(e0, e1));
        }
    }
    performanceLog << timestamp << "Finding edge journey pairs ends." << endl;
    // cout << "Found " << edgeJourneyPairs.size() << " edge journey pairs." << endl;

    // Count how many time each pair of edges appear.
    performanceLog << timestamp << "Deduplicate edge journey pairs begins." << endl;
    vector<uint64_t> coverage;
    deduplicateAndCountWithThreshold(edgeJourneyPairs, coverage, minTransitionGraphEdgeCoverage);
    edgeJourneyPairs.shrink_to_fit();
    // cout << "After  deduplicateAndCountWithThreshold there are " << edgeJourneyPairs.size() << " edge journey pairs." << endl;
    performanceLog << timestamp << "Deduplicate edge journey pairs ends." << endl;

    // Each of these generates an edge.
    for(uint64_t i=0; i<edgeJourneyPairs.size(); i++) {
        const auto& p = edgeJourneyPairs[i];
        const vertex_descriptor v0 = vertexMap[p.first];
        const vertex_descriptor v1 = vertexMap[p.second];
        add_edge(v0, v1, TransitionGraphEdge(coverage[i]), transitionGraph);
    }

    cout << "The TransitionGraph has " << num_vertices(transitionGraph) <<
        " vertices and " << num_edges(transitionGraph) << " edges." << endl;
}



void TransitionGraph::writeGraphviz(const string& fileName) const
{
    ofstream file(fileName);
    writeGraphviz(file);
}



void TransitionGraph::writeGraphviz(ostream& dot) const
{
    const TransitionGraph& transitionGraph = *this;

    dot << "digraph TransitionGraph {\n";

    BGL_FORALL_VERTICES(v, transitionGraph, TransitionGraph) {
        dot << transitionGraph[v].vertexId << ";\n";
    }

    BGL_FORALL_EDGES(e, transitionGraph, TransitionGraph) {
        const vertex_descriptor v0 = source(e, transitionGraph);
        const vertex_descriptor v1 = target(e, transitionGraph);

        dot << transitionGraph[v0].vertexId << "->";
        dot << transitionGraph[v1].vertexId << ";\n";
    }

    dot << "}\n";
}
