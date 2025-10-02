// Read following in the AssemblyGraph.

// Shasta.
#include "ReadFollowing.hpp"
#include "deduplicate.hpp"
#include "Journeys.hpp"
using namespace shasta;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>

// Standard library.
#include "fstream.hpp"
#include <iomanip>



ReadFollowing::ReadFollowing(
    const AssemblyGraph& assemblyGraph,
    uint64_t representativeRegionLength) :
    assemblyGraph(assemblyGraph)
{
    findAppearances(representativeRegionLength);
    countAppearances();
    findEdgePairs();
    createEdgePairsGraph();
    writeEdgePairsGraph();
}



#if 0
void ReadFollowing::createLineGraph()
{
    // Create a vertex for each AssemblyGraph edge.
    BGL_FORALL_EDGES(ae, assemblyGraph, AssemblyGraph) {
        lineGraphVertexMap.insert(make_pair(ae, add_vertex(ae, lineGraph)));
    }

    // Create edges of the LineGraph.
    BGL_FORALL_EDGES(ae0, assemblyGraph, AssemblyGraph) {
        const LineGraph::vertex_descriptor lv0 = lineGraphVertexMap[ae0];
        const AVertex av0 = target(ae0, assemblyGraph);
        BGL_FORALL_OUTEDGES(av0, ae1, assemblyGraph, AssemblyGraph) {
            const LineGraph::vertex_descriptor lv1 = lineGraphVertexMap[ae1];
            add_edge(lv0, lv1, lineGraph);
        }
    }

#if 0
    cout << "The AssemblyGraph has " << num_vertices(assemblyGraph) <<
        " vertices and " << num_edges(assemblyGraph) << " edges." << endl;
    cout << "The LineGraph has " << num_vertices(lineGraph) <<
        " vertices and " << num_edges(lineGraph) << " edges." << endl;
#endif
}
#endif



void ReadFollowing::findAppearances(uint64_t representativeRegionLength)
{
    const uint64_t orientedReadCount = assemblyGraph.journeys.size();
    initialAppearances.resize(orientedReadCount);
    finalAppearances.resize(orientedReadCount);



    // Loop over all edges of the AssemblyGraph.
    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        const AssemblyGraphEdge& edge = assemblyGraph[e];
        const uint64_t stepCount = edge.size();

        // Locate the initial representative region.
        const uint64_t initialBegin = 0;
        const uint64_t initialEnd = min(stepCount, representativeRegionLength);

        // Locate the final representative region.
        const uint64_t finalEnd = stepCount;
        const uint64_t finalBegin =
            ((stepCount >= representativeRegionLength) ? (stepCount - representativeRegionLength) : 0);

        // Appearances in the initial representative region of this edge.
        // For each OrientedReadId. store the last appearance in journey order.
        std::map<OrientedReadId, vector<uint32_t> > initialAppearancesMap;
        for(uint64_t stepId=initialBegin; stepId!=initialEnd; stepId++) {
            const AssemblyGraphEdgeStep& step = edge[stepId];
            const AnchorId anchorId = step.anchorPair.anchorIdA;
            for(const OrientedReadId orientedReadId: step.anchorPair.orientedReadIds) {
                const uint32_t positionInJourney =
                    assemblyGraph.anchors.getPositionInJourney(anchorId, orientedReadId);
                initialAppearancesMap[orientedReadId].push_back(positionInJourney);
            }
        }
        for(auto& p: initialAppearancesMap) {
            const OrientedReadId orientedReadId = p.first;
            vector<uint32_t>& positionsInJourney = p.second;
            std::ranges::sort(positionsInJourney);
            initialAppearances[orientedReadId.getValue()].push_back(
                Appearance(e, positionsInJourney.back()));
        }

        // Appearances in the final representative region of this edge.
        // For each OrientedReadId. store the first appearance in journey order.
        std::map<OrientedReadId, vector<uint32_t> > finalAppearancesMap;
        for(uint64_t stepId=finalBegin; stepId!=finalEnd; stepId++) {
            const AssemblyGraphEdgeStep& step = edge[stepId];
            const AnchorId anchorId = step.anchorPair.anchorIdB;
            for(const OrientedReadId orientedReadId: step.anchorPair.orientedReadIds) {
                const uint32_t positionInJourney =
                    assemblyGraph.anchors.getPositionInJourney(anchorId, orientedReadId);
                finalAppearancesMap[orientedReadId].push_back(positionInJourney);
            }
        }
        for(auto& p: finalAppearancesMap) {
            const OrientedReadId orientedReadId = p.first;
            vector<uint32_t>& positionsInJourney = p.second;
            std::ranges::sort(positionsInJourney);
            finalAppearances[orientedReadId.getValue()].push_back(
                Appearance(e, positionsInJourney.front()));
        }
    }



    // For each OrientedReadId, sort the appearances by edge id.
    class AppearanceSorter {
    public:
        const AssemblyGraph& assemblyGraph;
        AppearanceSorter(const AssemblyGraph& assemblyGraph) : assemblyGraph(assemblyGraph) {}
        bool operator()(const Appearance& x, const Appearance& y) const
        {
            return assemblyGraph[x.e].id < assemblyGraph[y.e].id;
        }

    };
    for(vector<Appearance>& v: initialAppearances) {
        sort(v.begin(), v.end(), AppearanceSorter(assemblyGraph));
    }
    for(vector<Appearance>& v: finalAppearances) {
        sort(v.begin(), v.end(), AppearanceSorter(assemblyGraph));
    }
}



void ReadFollowing::countAppearances()
{
    const uint64_t orientedReadCount = assemblyGraph.journeys.size();


    // Loop over all OrientedReadIds.
    for(uint64_t i=0; i<orientedReadCount; i++) {

        // Initial appearances.
        for(const Appearance appearance: initialAppearances[i]) {
            const AEdge e = appearance.e;
            const auto it = initialAppearancesCount.find(e);
            if(it == initialAppearancesCount.end()) {
                initialAppearancesCount.insert(make_pair(e, 1));
            } else {
                ++(it->second);
            }
        }

        // Final appearances.
        for(const Appearance appearance: finalAppearances[i]) {
            const AEdge e = appearance.e;
            const auto it = finalAppearancesCount.find(e);
            if(it == finalAppearancesCount.end()) {
                finalAppearancesCount.insert(make_pair(e, 1));
            } else {
                ++(it->second);
            }
        }

    }

}



void ReadFollowing::findEdgePairs()
{
    const uint64_t orientedReadCount = assemblyGraph.journeys.size();

    vector<AssemblyGraphEdgePair> edgePairsVector;

    // Loop over all OrientedReadIds.
    for(uint64_t i=0; i<orientedReadCount; i++) {

        // Look over pairs (final appearance, initial appearance).
        for(const Appearance appearance0: finalAppearances[i]) {
            const AEdge e0 = appearance0.e;
            for(const Appearance appearance1: initialAppearances[i]) {
                const AEdge e1 = appearance1.e;
                if((e1 != e0) and (appearance0.positionInJourney < appearance1.positionInJourney)) {
                    edgePairsVector.push_back(make_pair(e0, e1));
                }
            }
        }
    }
    vector<uint64_t> coverage;
    deduplicateAndCount(edgePairsVector, coverage);
    SHASTA_ASSERT(coverage.size() == edgePairsVector.size());

    // Now store them in the edgePairs map.
    for(uint64_t i=0; i<edgePairsVector.size(); i++) {
        edgePairs.insert(make_pair(edgePairsVector[i], coverage[i]));
    }

    // Write out the edgePairs.
    {
        ofstream csv("ReadFollowing-EdgePairs.csv");
        csv << "Id0,Id1,Coverage,\n";
        for(const auto& p: edgePairs) {
            const auto& edgePair = p.first;
            const uint64_t coverage = p.second;
            const AEdge e0 = edgePair.first;
            const AEdge e1 = edgePair.second;
            csv << assemblyGraph[e0].id << ",";
            csv << assemblyGraph[e1].id << ",";
            csv << coverage << ",";
            csv << "\n";
        }
    }
}



#if 0
void ReadFollowing::createEdgePairsGraph()
{
    // Create a vertex for each AssemblyGraph edge.
    BGL_FORALL_EDGES(ae, assemblyGraph, AssemblyGraph) {
        edgePairsVertexMap.insert(make_pair(ae, add_vertex(ae, edgePairsGraph)));
    }

    // Each AssemblyGraphEdgePair with sufficient coverage generates an edge.
    for(const auto& p: edgePairs) {
        const AssemblyGraphEdgePair& edgePair = p.first;
        const uint64_t coverage = p.second;
        if(coverage >= minCoverage) {
            const AEdge ae0 = edgePair.first;
            const AEdge ae1 = edgePair.second;
            const EdgePairsGraph::vertex_descriptor v0 = edgePairsVertexMap[ae0];
            const EdgePairsGraph::vertex_descriptor v1 = edgePairsVertexMap[ae1];
            add_edge(v0, v1, coverage, edgePairsGraph);
        }
    }
}
#endif



void ReadFollowing::writeEdgePairsGraph()
{
    ofstream dot("ReadFollowing.dot");
    dot << "digraph EdgePairsGraph {\n";
    dot << std::fixed << std::setprecision(1);

    BGL_FORALL_VERTICES(v, edgePairsGraph, EdgePairsGraph) {
        const EdgePairsGraphVertex& edgePairsGraphVertex = edgePairsGraph[v];
        const AEdge ae = edgePairsGraphVertex.ae;
        const AssemblyGraphEdge& assemblyGraphEdge = assemblyGraph[ae];
        dot << assemblyGraphEdge.id <<
            " [label=\"" << assemblyGraphEdge.id << "\\n" <<
            edgePairsGraphVertex.length << "\\n" <<
            getInitialAppearancesCount(ae) << "/" <<
            getFinalAppearancesCount(ae) <<
            "\"]"
            ";\n";
    }

    BGL_FORALL_EDGES(e, edgePairsGraph, EdgePairsGraph) {
        const EdgePairsGraph::vertex_descriptor v0 = source(e, edgePairsGraph);
        const EdgePairsGraph::vertex_descriptor v1 = target(e, edgePairsGraph);
        const AEdge ae0 = edgePairsGraph[v0].ae;
        const AEdge ae1 = edgePairsGraph[v1].ae;
        const uint64_t coverage = edgePairsGraph[e];
        dot <<
            assemblyGraph[ae0].id << "->" <<
            assemblyGraph[ae1].id <<
            " [penwidth=\"" << 0.4 * double(coverage) << "\""
            " tooltip=\"" << coverage << "\"]"
            ";\n";
    }

    dot << "}\n";
}



#if 0
void ReadFollowing::followForward(AEdge ae0) const
{
    // EXPOSE WHEN CODE STABILIZES.
    const uint64_t minCoverage = 3;
    const double minCoverageFraction = 0.8;

    const auto it0 = edgePairsVertexMap.find(ae0);
    SHASTA_ASSERT(it0 != edgePairsVertexMap.end());
    const EdgePairsGraph::vertex_descriptor v0 = it0->second;

    const auto it = finalAppearancesCount.find(ae0);
    cout << assemblyGraph[ae0].id << " " << it->second << endl;

    // Find possible following AEdges, grouping them by they source AVertex.
    std::map<AVertex, vector< pair<AEdge, uint64_t> > > m;
    BGL_FORALL_OUTEDGES(v0, e, edgePairsGraph, EdgePairsGraph) {
        const EdgePairsGraph::vertex_descriptor v1 = target(e, edgePairsGraph);
        const AEdge ae1 = edgePairsGraph[v1];
        const uint64_t coverage = edgePairsGraph[e];
        const AVertex av1 = source(ae1, assemblyGraph);
        m[av1].push_back(make_pair(ae1, coverage));
    }


    // Loop over the source vertices.
    std::set<AEdge> next;
    for(const auto& p: m) {
        const AVertex v1 = p.first;
        const vector< pair<AEdge, uint64_t> >& v = p.second;

        // Compute total coverage.
        uint64_t totalCoverage = 0;
        for(const auto& q: v) {
            const uint64_t coverage = q.second;
            totalCoverage += coverage;
        }

        cout << assemblyGraph[v1].id << ":";
        for(const auto& q: v) {
            const AEdge ae = q.first;
            const uint64_t coverage = q.second;
            cout << " (" << assemblyGraph[ae].id << "," << coverage;
            if((coverage >= minCoverage) and ((double(coverage) / double(totalCoverage)) >= minCoverageFraction)) {
                cout << " ***";
                next.insert(ae);
            }
            cout << ")";
        }
        cout << endl;
    }

    ofstream csv("ReadFollowing-Bandage.csv");
    csv << "Segment,Color\n";
    BGL_FORALL_EDGES(ae, assemblyGraph, AssemblyGraph) {
        string color;
        if(ae == ae0) {
            color = "Blue";
        } else {
            if(next.contains(ae)) {
                color = "Green";
            } else {
                color = "LightGrey";
            }
        }
        csv << assemblyGraph[ae].id << "," << color << "\n";
    }
}



void ReadFollowing::followBackward(AEdge ae0) const
{
    // EXPOSE WHEN CODE STABILIZES.
    const uint64_t minCoverage = 3;
    const double minCoverageFraction = 0.8;

    const auto it0 = edgePairsVertexMap.find(ae0);
    SHASTA_ASSERT(it0 != edgePairsVertexMap.end());
    const EdgePairsGraph::vertex_descriptor v0 = it0->second;

    const auto it = finalAppearancesCount.find(ae0);
    cout << assemblyGraph[ae0].id << " " << it->second << endl;

    // Find possible previous AEdges, grouping them by they target AVertex.
    std::map<AVertex, vector< pair<AEdge, uint64_t> > > m;
    BGL_FORALL_INEDGES(v0, e, edgePairsGraph, EdgePairsGraph) {
        const EdgePairsGraph::vertex_descriptor v1 = source(e, edgePairsGraph);
        const AEdge ae1 = edgePairsGraph[v1];
        const uint64_t coverage = edgePairsGraph[e];
        const AVertex av1 = target(ae1, assemblyGraph);
        m[av1].push_back(make_pair(ae1, coverage));
    }


    // Loop over the target vertices.
    std::set<AEdge> next;
    for(const auto& p: m) {
        const AVertex v1 = p.first;
        const vector< pair<AEdge, uint64_t> >& v = p.second;

        // Compute total coverage.
        uint64_t totalCoverage = 0;
        for(const auto& q: v) {
            const uint64_t coverage = q.second;
            totalCoverage += coverage;
        }

        cout << assemblyGraph[v1].id << ":";
        for(const auto& q: v) {
            const AEdge ae = q.first;
            const uint64_t coverage = q.second;
            cout << " (" << assemblyGraph[ae].id << "," << coverage;
            if((coverage >= minCoverage) and ((double(coverage) / double(totalCoverage)) >= minCoverageFraction)) {
                cout << " ***";
                next.insert(ae);
            }
            cout << ")";
        }
        cout << endl;
    }

    ofstream csv("ReadFollowing-Bandage.csv");
    csv << "Segment,Color\n";
    BGL_FORALL_EDGES(ae, assemblyGraph, AssemblyGraph) {
        string color;
        if(ae == ae0) {
            color = "Blue";
        } else {
            if(next.contains(ae)) {
                color = "Green";
            } else {
                color = "LightGrey";
            }
        }
        csv << assemblyGraph[ae].id << "," << color << "\n";
    }
}



void ReadFollowing::followForward(
    AEdge ae0,
    vector<AssemblyGraph::edge_descriptor>& path) const
{
    const auto it0 = edgePairsVertexMap.find(ae0);
    SHASTA_ASSERT(it0 != edgePairsVertexMap.end());
    const EdgePairsGraph::vertex_descriptor v0 = it0->second;

    std::set<AEdge> next;
    next.insert(ae0);
    BGL_FORALL_OUTEDGES(v0, e, edgePairsGraph, EdgePairsGraph) {
        const EdgePairsGraph::vertex_descriptor v1 = target(e, edgePairsGraph);
        const AEdge e1 = edgePairsGraph[v1];
        next.insert(e1);
    }

    ofstream csv("ReadFollowing-Bandage.csv");
    csv << "Segment,Color\n";
    BGL_FORALL_EDGES(ae, assemblyGraph, AssemblyGraph) {
        string color;
        if(ae == ae0) {
            color = "Blue";
        } else {
            if(next.contains(ae)) {
                color = "Green";
            } else {
                color = "LightGrey";
            }
        }
        csv << assemblyGraph[ae].id << "," << color << "\n";
    }



    const FilteringPredicate predicate(&lineGraph, &next);
    FilteredLineGraph filteredLineGraph(lineGraph, predicate, predicate);

    // Write out the FilteredLineGraph.
    ofstream dot("ReadFollowing-FilteredLineGraph.dot");
    dot << "digraph FilteredLineGraph {\n";
    for(const AEdge ae: next) {
        dot << assemblyGraph[ae].id;
        uint64_t coverage;
        if(ae == ae0) {
            const auto it = finalAppearancesCount.find(ae);
            SHASTA_ASSERT(it != finalAppearancesCount.end());
            coverage = it->second;
        } else {
            const auto it = edgePairs.find(AssemblyGraphEdgePair(ae0, ae));
            SHASTA_ASSERT(it != edgePairs.end());
            coverage = it->second;
        }
        dot << " [label=\"" << assemblyGraph[ae].id << "\\n" << coverage << "\"]";
        dot << ";\n";
    }
    for(const AEdge ae0: next) {
        const auto it0 = lineGraphVertexMap.find(ae0);
        SHASTA_ASSERT(it0 != lineGraphVertexMap.end());
        const LineGraph::vertex_descriptor lv0 = it0->second;
        BGL_FORALL_OUTEDGES(lv0, le, filteredLineGraph, FilteredLineGraph) {
            const LineGraph::vertex_descriptor lv1 = target(le, filteredLineGraph);
            const AEdge ae1 = filteredLineGraph[lv1];
            dot << assemblyGraph[ae0].id << "->" <<
                assemblyGraph[ae1].id << ";\n";
        }
    }
    dot << "}\n";



    // Follow the FilteredLineGraph forward starting at the vertex corresponding to ae0.
    // Stop when the out-degree is not exactly one.
    path.clear();
    const auto it = lineGraphVertexMap.find(ae0);
    SHASTA_ASSERT(it != lineGraphVertexMap.end());
    FilteredLineGraph::vertex_descriptor lv = it->second;
    while(true) {
        const AEdge ae = lineGraph[lv];
        path.push_back(ae);
        if(out_degree(lv, filteredLineGraph) != 1) {
            break;
        }
        FilteredLineGraph::out_edge_iterator it;
        tie(it, ignore) = out_edges(lv, filteredLineGraph);
        lv = target(*it, filteredLineGraph);
    }

    for(const AEdge& ae: path) {
        cout << assemblyGraph[ae].id << " ";
    }
    cout << endl;

}



void ReadFollowing::followBackward(
    AEdge ae0,
    vector<AssemblyGraph::edge_descriptor>& path) const
{
    const auto it0 = edgePairsVertexMap.find(ae0);
    SHASTA_ASSERT(it0 != edgePairsVertexMap.end());
    const EdgePairsGraph::vertex_descriptor v0 = it0->second;

    std::set<AEdge> next;
    next.insert(ae0);
    BGL_FORALL_INEDGES(v0, e, edgePairsGraph, EdgePairsGraph) {
        const EdgePairsGraph::vertex_descriptor v1 = source(e, edgePairsGraph);
        const AEdge e1 = edgePairsGraph[v1];
        next.insert(e1);
    }

    ofstream csv("ReadFollowing-Bandage.csv");
    csv << "Segment,Color\n";
    BGL_FORALL_EDGES(ae, assemblyGraph, AssemblyGraph) {
        string color;
        if(ae == ae0) {
            color = "Blue";
        } else {
            if(next.contains(ae)) {
                color = "Green";
            } else {
                color = "LightGrey";
            }
        }
        csv << assemblyGraph[ae].id << "," << color << "\n";
    }



    const FilteringPredicate predicate(&lineGraph, &next);
    FilteredLineGraph filteredLineGraph(lineGraph, predicate, predicate);

    // Write out the FilteredLineGraph.
    ofstream dot("ReadFollowing-FilteredLineGraph.dot");
    dot << "digraph FilteredLineGraph {\n";
    for(const AEdge ae: next) {
        dot << assemblyGraph[ae].id;
        uint64_t coverage;
        if(ae == ae0) {
            const auto it = initialAppearancesCount.find(ae);
            SHASTA_ASSERT(it != initialAppearancesCount.end());
            coverage = it->second;
        } else {
            const auto it = edgePairs.find(AssemblyGraphEdgePair(ae, ae0));
            SHASTA_ASSERT(it != edgePairs.end());
            coverage = it->second;
        }
        dot << " [label=\"" << assemblyGraph[ae].id << "\\n" << coverage << "\"]";
        dot << ";\n";
    }
    for(const AEdge ae0: next) {
        const auto it0 = lineGraphVertexMap.find(ae0);
        SHASTA_ASSERT(it0 != lineGraphVertexMap.end());
        const LineGraph::vertex_descriptor lv0 = it0->second;
        BGL_FORALL_OUTEDGES(lv0, le, filteredLineGraph, FilteredLineGraph) {
            const LineGraph::vertex_descriptor lv1 = target(le, filteredLineGraph);
            const AEdge ae1 = filteredLineGraph[lv1];
            dot << assemblyGraph[ae0].id << "->" <<
                assemblyGraph[ae1].id << ";\n";
        }
    }
    dot << "}\n";



    // Follow the FilteredLineGraph backward starting at the vertex corresponding to ae0.
    // Stop when the out-degree is not exactly one.
    path.clear();
    const auto it = lineGraphVertexMap.find(ae0);
    SHASTA_ASSERT(it != lineGraphVertexMap.end());
    FilteredLineGraph::vertex_descriptor lv = it->second;
    while(true) {
        const AEdge ae = lineGraph[lv];
        path.push_back(ae);
        if(in_degree(lv, filteredLineGraph) != 1) {
            break;
        }
        FilteredLineGraph::in_edge_iterator it;
        tie(it, ignore) = in_edges(lv, filteredLineGraph);
        lv = source(*it, filteredLineGraph);
    }
    std::ranges::reverse(path);

    for(const AEdge& ae: path) {
        cout << assemblyGraph[ae].id << " ";
    }
    cout << endl;
}
#endif



void ReadFollowing::createEdgePairsGraph()
{
    // EXPOSE WHEN CODE STABILIZES.
    const uint64_t minCoverage = 3;
    const double minCoverageFraction = 0.8;
    const uint64_t maxAppearanceCount = 25;

    // Create a vertex for each AssemblyGraph edge.
    BGL_FORALL_EDGES(ae, assemblyGraph, AssemblyGraph) {
        edgePairsVertexMap.insert(make_pair(ae,
            add_vertex(EdgePairsGraphVertex(assemblyGraph, ae), edgePairsGraph)));
    }

    // Group AssemblyGraphEdgePairs (e0, e1)
    // by (v0, v1), where v0 = target(e0), v1 = source(e1).
    std::set<AVertexPair> s;
    for(const auto& p: edgePairs) {
        const AEdgePair& aEdgePair = p.first;
        const AEdge e0 = aEdgePair.first;
        const AEdge e1 = aEdgePair.second;
        const AVertex v0 = target(e0, assemblyGraph);
        const AVertex v1 = source(e1, assemblyGraph);
        s.insert(make_pair(v0, v1));
    }

    // Loop over pairs of AssemblyGraph vertices (v0, v1) such that there
    // is at least one EdgePair (e0, e1) with v0 = target(e0, assemblyGraph)
    // and v1 = source(e1, assemblyGraph).
    for(const auto& p: s) {
        const AVertex v0 = p.first;
        const AVertex v1 = p.second;

        // Gather the in-edges of v0.
        vector<AEdge> in;
        BGL_FORALL_INEDGES(v0, e, assemblyGraph, AssemblyGraph) {
            in.push_back(e);
        }
        sort(in.begin(), in.end(), assemblyGraph.orderById);
        const uint64_t nIn = in.size();

        // Gather the out-edges of v1.
        vector<AEdge> out;
        BGL_FORALL_OUTEDGES(v1, e, assemblyGraph, AssemblyGraph) {
            out.push_back(e);
        }
        sort(out.begin(), out.end(), assemblyGraph.orderById);
        const uint64_t nOut = out.size();

#if 0
        cout << "In:";
        for(const AEdge e: in) {
            cout << " " << assemblyGraph[e].id;
        }
        cout << endl;

        cout << "Out:";
        for(const AEdge e: out) {
            cout << " " << assemblyGraph[e].id;
        }
        cout << endl;
#endif
        // Create a tangle matrix using the edgePairs.
        vector< vector<uint64_t> > tangleMatrix(nIn, vector<uint64_t>(nOut));
        // cout << "Tangle matrix:" << endl;
        for(uint64_t i0=0; i0<nIn; i0++) {
            const AEdge e0 = in[i0];
            for(uint64_t i1=0; i1<nOut; i1++) {
                const AEdge e1 = out[i1];
                uint64_t coverage = 0;
                auto it = edgePairs.find(make_pair(e0, e1));
                if(it != edgePairs.end()) {
                    coverage = it->second;
                }
                // cout << coverage << ",";
                tangleMatrix[i0][i1] = coverage;
            }
            // cout << endl;
        }

        // Compute the sum of tangle matrix values for each of the in-edges.
        vector<uint64_t> inSum(nIn, 0);
        for(uint64_t i0=0; i0<nIn; i0++) {
            for(uint64_t i1=0; i1<nOut; i1++) {
                inSum[i0] += tangleMatrix[i0][i1];
            }
        }
        /*
        cout << "inSum: ";
        for(const uint64_t s: inSum) {
            cout << s << ",";
        }
        cout << endl;
        */

        // Compute the sum of tangle matrix values for each of the out-edges.
        vector<uint64_t> outSum(nOut, 0);
        for(uint64_t i0=0; i0<nIn; i0++) {
            for(uint64_t i1=0; i1<nOut; i1++) {
                outSum[i1] += tangleMatrix[i0][i1];
            }
        }
        /*
        cout << "outSum: ";
        for(const uint64_t s: outSum) {
            cout << s << ",";
        }
        cout << endl;
        */



        // Each "strong" element of the tangle matrix generates a strong edge pair
        // and an edge of the EdgePairsGraph.
        // "Strong" means:
        // - coverage >= minCoverage
        // - coverage / sum(tangle matrix values for same in-edge) >= minCoverageFraction.
        // - coverage / sum(tangle matrix values for same out-edge) >= minCoverageFraction.
        // - The number of final/initial appearances is no more than maxAppearanceCount.
        for(uint64_t i0=0; i0<nIn; i0++) {
            for(uint64_t i1=0; i1<nOut; i1++) {
                const uint64_t coverage = tangleMatrix[i0][i1];
                if(coverage < minCoverage) {
                    continue;
                }
                if(double(coverage) / double(inSum[i0]) < minCoverageFraction) {
                    continue;
                }
                if(double(coverage) / double(outSum[i1]) < minCoverageFraction) {
                    continue;
                }
                const AEdge e0 = in[i0];
                const AEdge e1 = out[i1];

                if(finalAppearancesCount[e0] > maxAppearanceCount) {
                    continue;
                }
                if(initialAppearancesCount[e1] > maxAppearanceCount) {
                    continue;
                }

                add_edge(edgePairsVertexMap[e0], edgePairsVertexMap[e1], coverage, edgePairsGraph);
            }
        }
    }
}



// In the EdgePairsGraph, find a path that starts at a given AEdge
// and moves forward/backward. At each step we choose the child vertex
// corresponding to the longest AEdge.
void ReadFollowing::findPath(AEdge ae, uint64_t direction) const
{
    if(direction == 0) {
        findForwardPath(ae);
    } else if(direction == 1) {
        findBackwardPath(ae);
    } else {
        SHASTA_ASSERT(0);
    }
}



void ReadFollowing::findForwardPath(AEdge ae) const
{
    const auto it = edgePairsVertexMap.find(ae);
    SHASTA_ASSERT(it != edgePairsVertexMap.end());
    EdgePairsGraph::vertex_descriptor v = it->second;

    vector<AEdge> path;
    while(true) {
        path.push_back(edgePairsGraph[v].ae);

        EdgePairsGraph::vertex_descriptor vNext = EdgePairsGraph::null_vertex();
        uint64_t bestLength = 0;
        BGL_FORALL_OUTEDGES(v, e, edgePairsGraph, EdgePairsGraph) {
            EdgePairsGraph::vertex_descriptor v1 = target(e, edgePairsGraph);
            const uint64_t length = edgePairsGraph[v1].length;
            if(length > bestLength) {
                vNext = v1;
                bestLength = length;
            }
        }

        if(vNext == EdgePairsGraph::null_vertex()) {
            break;
        }

        v = vNext;
    }

    cout << "Path:" << endl;
    for(const AEdge ae: path) {
        cout << assemblyGraph[ae].id << ",";
    }
    cout << endl;
}



void ReadFollowing::findBackwardPath(AEdge ae) const
{
    const auto it = edgePairsVertexMap.find(ae);
    SHASTA_ASSERT(it != edgePairsVertexMap.end());
    EdgePairsGraph::vertex_descriptor v = it->second;

    vector<AEdge> path;
    while(true) {
        path.push_back(edgePairsGraph[v].ae);

        EdgePairsGraph::vertex_descriptor vNext = EdgePairsGraph::null_vertex();
        uint64_t bestLength = 0;
        BGL_FORALL_INEDGES(v, e, edgePairsGraph, EdgePairsGraph) {
            EdgePairsGraph::vertex_descriptor v1 = source(e, edgePairsGraph);
            const uint64_t length = edgePairsGraph[v1].length;
            if(length > bestLength) {
                vNext = v1;
                bestLength = length;
            }
        }

        if(vNext == EdgePairsGraph::null_vertex()) {
            break;
        }

        v = vNext;
    }
    std::ranges::reverse(path);

    cout << "Path:" << endl;
    for(const AEdge ae: path) {
        cout << assemblyGraph[ae].id << ",";
    }
    cout << endl;
}



ReadFollowing::EdgePairsGraphVertex::EdgePairsGraphVertex(
    const AssemblyGraph& assemblyGraph,
    AEdge ae) :
    ae(ae)
{
    const AssemblyGraphEdge& edge = assemblyGraph[ae];
    if(edge.wasAssembled) {
        length = edge.sequenceLength();
    } else {
        length = edge.offset();
    }
}



uint64_t ReadFollowing::getInitialAppearancesCount(AEdge ae) const
{
    const auto it = initialAppearancesCount.find(ae);
    if(it == initialAppearancesCount.end()) {
        return 0;
    } else {
        return it->second;
    }
}



uint64_t ReadFollowing::getFinalAppearancesCount(AEdge ae) const
{
    const auto it = finalAppearancesCount.find(ae);
    if(it == finalAppearancesCount.end()) {
        return 0;
    } else {
        return it->second;
    }
}

