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



ReadFollowing::ReadFollowing(const AssemblyGraph& assemblyGraph) :
    assemblyGraph(assemblyGraph)
{
    findAppearances();
    countAppearances();
    findEdgePairs();
    createEdgePairsGraph();
    writeEdgePairsGraph();
}



void ReadFollowing::findAppearances()
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



void ReadFollowing::writeEdgePairsGraph()
{
    ofstream dot("ReadFollowing.dot");
    dot << "digraph EdgePairsGraph {\n";
    dot << std::fixed << std::setprecision(2);

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

        const double j = jaccard(e);
        SHASTA_ASSERT(j <= 1.);
        const double hue = j / 3.; // So 0=red, 1=green.

        dot <<
            assemblyGraph[ae0].id << "->" <<
            assemblyGraph[ae1].id <<
            " [penwidth=\"" << 0.4 * double(coverage) << "\""
            " tooltip=\"" << coverage << " " << j << "\""
            " color=\"" << hue << ",1,.9\""
            "]"
            ";\n";
    }

    dot << "}\n";
}



void ReadFollowing::createEdgePairsGraph()
{

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

                if(initialAppearancesCount[e0] > maxAppearanceCount) {
                    continue;
                }
                if(finalAppearancesCount[e0] > maxAppearanceCount) {
                    continue;
                }
                if(initialAppearancesCount[e1] > maxAppearanceCount) {
                    continue;
                }
                if(finalAppearancesCount[e1] > maxAppearanceCount) {
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
        // uint64_t bestLength = 0;
        double bestJaccard = 0.;
        BGL_FORALL_OUTEDGES(v, e, edgePairsGraph, EdgePairsGraph) {
            EdgePairsGraph::vertex_descriptor v1 = target(e, edgePairsGraph);
            // const uint64_t length = edgePairsGraph[v1].length;
            const double j = jaccard(e);
            if(j > bestJaccard /*length > bestLength */) {
                vNext = v1;
                // bestLength = length;
                bestJaccard = j;
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
        // uint64_t bestLength = 0;
        double bestJaccard = 0.;
        BGL_FORALL_INEDGES(v, e, edgePairsGraph, EdgePairsGraph) {
            EdgePairsGraph::vertex_descriptor v1 = source(e, edgePairsGraph);
            // const uint64_t length = edgePairsGraph[v1].length;
            const double j = jaccard(e);
            if(j > bestJaccard /*length > bestLength */) {
                vNext = v1;
                // bestLength = length;
                bestJaccard = j;
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



// Jaccard similarity for an EdgePairsGraph.
// Computed using the finalAppearancesCount of the surce vertex
// and the initialAppearancesCount of the target vertex.
double ReadFollowing::jaccard(EdgePairsGraph::edge_descriptor e) const
{
    const EdgePairsGraph::vertex_descriptor v0 = source(e, edgePairsGraph);
    const EdgePairsGraph::vertex_descriptor v1 = target(e, edgePairsGraph);

    const AEdge ae0 = edgePairsGraph[v0].ae;
    const AEdge ae1 = edgePairsGraph[v1].ae;

    const uint64_t n0 = getFinalAppearancesCount(ae0);
    const uint64_t n1 = getInitialAppearancesCount(ae1);

    const uint64_t n01 = edgePairsGraph[e];

    const uint64_t intersectionSize = n01;
    const uint64_t unionSize = n0 + n1 - intersectionSize;

    return double(intersectionSize) / double(unionSize);
}
