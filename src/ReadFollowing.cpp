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



ReadFollowing::ReadFollowing(
    const AssemblyGraph& assemblyGraph,
    uint64_t representativeRegionLength,
    uint64_t minCoverage) :
    assemblyGraph(assemblyGraph)
{
    findAppearances(representativeRegionLength);
    findEdgePairs();
    createEdgePairsGraph(minCoverage);
    writeEdgePairsGraph();
}



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



void ReadFollowing::createEdgePairsGraph(uint64_t minCoverage)
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



void ReadFollowing::writeEdgePairsGraph()
{
    ofstream dot("ReadFollowing-EdgePairsGraph.dot");
    dot << "digraph EdgePairsGraph {\n";

    BGL_FORALL_VERTICES(v, edgePairsGraph, EdgePairsGraph) {
        const AEdge ae = edgePairsGraph[v];
        dot << assemblyGraph[ae].id << ";\n";
    }

    BGL_FORALL_EDGES(e, edgePairsGraph, EdgePairsGraph) {
        const EdgePairsGraph::vertex_descriptor v0 = source(e, edgePairsGraph);
        const EdgePairsGraph::vertex_descriptor v1 = target(e, edgePairsGraph);
        const AEdge ae0 = edgePairsGraph[v0];
        const AEdge ae1 = edgePairsGraph[v1];
        const uint64_t coverage = edgePairsGraph[e];
        dot <<
            assemblyGraph[ae0].id << "->" <<
            assemblyGraph[ae1].id <<
            " [label=\"" << coverage << "\"]"
            ";\n";
    }

    dot << "}\n";
}



void ReadFollowing::followForward(AEdge ae0) const
{
    const auto it0 = edgePairsVertexMap.find(ae0);
    SHASTA_ASSERT(it0 != edgePairsVertexMap.end());
    const EdgePairsGraph::vertex_descriptor v0 = it0->second;

    std::set<AEdge> next;
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
}
