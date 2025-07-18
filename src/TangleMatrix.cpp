#include "TangleMatrix.hpp"
#include "Anchor.hpp"
#include "Base.hpp"
using namespace shasta;



TangleMatrix::TangleMatrix(
    const AssemblyGraph& assemblyGraph,
    vector<edge_descriptor> entranceEdges,
    vector<edge_descriptor> exitEdges,
    uint64_t maxTrim,
    double aDrift,
    double bDrift)
{
    const bool debug = false;

    SHASTA_ASSERT(std::ranges::is_sorted(entranceEdges, AssemblyGraph::OrderById(assemblyGraph)));
    for(const edge_descriptor e: entranceEdges) {
        const AssemblyGraphEdge& edge = assemblyGraph[e];

        // Among the last maxTrim steps, pick the one with the best coverage.
        uint64_t bestStepId = invalid<uint64_t>;
        uint64_t bestTrim = invalid<uint64_t>;
        uint64_t bestCoverage = 0;
        for(uint64_t trim=0; trim<=min(maxTrim, edge.size() - 1); trim++) {
            const uint64_t stepId = edge.size() - 1 - trim;
            const uint64_t coverage = edge[stepId].anchorPair.size();
            if(coverage > bestCoverage) {
                bestStepId = stepId;
                bestTrim = trim;
                bestCoverage = coverage;
            }
        }

        entrances.emplace_back(e, bestTrim, edge[bestStepId]);
        if(debug) {
            cout << "Entrance " << assemblyGraph[e].id << ": chose step " << bestStepId << " coverage " << bestCoverage << endl;
        }
    }



    SHASTA_ASSERT(std::ranges::is_sorted(exitEdges, AssemblyGraph::OrderById(assemblyGraph)));
    for(const edge_descriptor e: exitEdges) {
        const AssemblyGraphEdge& edge = assemblyGraph[e];

        // Among the first maxTrim steps, pick the one with the best coverage.
        uint64_t bestStepId = invalid<uint64_t>;
        uint64_t bestTrim = invalid<uint64_t>;
        uint64_t bestCoverage = 0;
        for(uint64_t trim=0; trim<=min(maxTrim, edge.size() - 1); trim++) {
            const uint64_t stepId = trim;
            const uint64_t coverage = edge[stepId].anchorPair.size();
            if(coverage > bestCoverage) {
                bestStepId = stepId;
                bestTrim = trim;
                bestCoverage = coverage;
            }
        }

        exits.emplace_back(e, bestTrim, edge[bestStepId]);
        if(debug) {
            cout << "Exit " << assemblyGraph[e].id << ": chose step " << bestStepId << " coverage " << bestCoverage << endl;
        }
    }



    // Compute the tangle matrix.
    tangleMatrix.resize(entrances.size(), vector<AnchorPair>(exits.size()));
    for(uint64_t iEntrance=0; iEntrance<entrances.size(); iEntrance++) {
        EntranceOrExit& entrance = entrances[iEntrance];
        for(uint64_t iExit=0; iExit<exits.size(); iExit++) {
            EntranceOrExit& exit = exits[iExit];

            tangleMatrix[iEntrance][iExit] = assemblyGraph.anchors.bridge(
                entrance.step.anchorPair,
                exit.step.anchorPair,
                aDrift, bDrift);
            const uint64_t coverage = tangleMatrix[iEntrance][iExit].orientedReadIds.size();
            entrance.commonCoverage += coverage;
            exit.commonCoverage += coverage;
        }
    }
}



void TangleMatrix::writeHtml(
    const AssemblyGraph& assemblyGraph,
    ostream& html) const
{
    html <<
        "<h2>Tangle matrix</h2>"
        "<p>"
        "<table>"
        "<tr><td rowspan=5 colspan=5>"
        "<th colspan=" << (exits.size() + 1) << " style='background-color:LightCyan'>Exits";

    html << "<tr><th style='background-color:LightCyan'>Index";
    for(uint64_t i=0; i<exits.size(); i++) {
        html << "<td class=centered style='background-color:LightCyan'>" << i;
    }

    html << "<tr><th style='background-color:LightCyan'>Segment";
    for(uint64_t i=0; i<exits.size(); i++) {
        html << "<td class=centered style='background-color:LightCyan'>" << assemblyGraph[exits[i].e].id;
    }

    html << "<tr><th style='background-color:LightCyan'>Coverage";
    for(uint64_t i=0; i<exits.size(); i++) {
        html << "<td class=centered style='background-color:LightCyan'>" << exits[i].coverage();
    }

    html << "<tr><th style='background-color:LightCyan'>Common<br>coverage";
    for(uint64_t i=0; i<exits.size(); i++) {
        html << "<td class=centered style='background-color:LightCyan'>" << exits[i].commonCoverage;
    }

    html <<
        "<tr><th rowspan=" << entrances.size() + 1 <<
        " style='background-color:CornSilk'>Entrances"
        "<th style='background-color:CornSilk'>Index"
        "<th style='background-color:CornSilk'>Segment"
        "<th style='background-color:CornSilk'>Coverage"
        "<th style='background-color:CornSilk'>Common<br>coverage<td>";

    for(uint64_t i=0; i<exits.size(); i++) {
        html << "<td>";
    }

    for(uint64_t iEntrance=0; iEntrance<entrances.size(); iEntrance++) {
        const EntranceOrExit& entrance = entrances[iEntrance];
        html <<
            "<tr>"
            "<td class=centered style='background-color:CornSilk'>" << iEntrance <<
            "<td class=centered style='background-color:CornSilk'>" << assemblyGraph[entrance.e].id <<
            "<td class=centered style='background-color:CornSilk'>" << entrance.coverage() <<
            "<td class=centered style='background-color:CornSilk'>" << entrance.commonCoverage <<
            "<td>";
        for(uint64_t iExit=0; iExit<exits.size(); iExit++) {
            html << "<td class=centered style='background-color:LightPink'>" <<
                tangleMatrix[iEntrance][iExit].orientedReadIds.size();
        }
    }

    html << "</table>";

}



// Read following on the entrances/exits.
void TangleMatrix::readFollowing(const AssemblyGraph& assemblyGraph)
{
    const bool debug = false;
    readFollowingMap.clear();

    for(uint64_t iEntrance=0; iEntrance<entrances.size(); iEntrance++) {
        const auto& entrance = entrances[iEntrance];
        const AssemblyGraph::edge_descriptor e = entrance.e;
        const AssemblyGraphEdge& edge = assemblyGraph[e];
        for(uint64_t stepId=0; stepId<edge.size(); stepId++) {
            const AssemblyGraphEdgeStep& step = edge[stepId];
            const AnchorPair& anchorPair = step.anchorPair;
            for(const OrientedReadId orientedReadId: anchorPair.orientedReadIds) {
                StepIdentifier stepIdentifier;
                stepIdentifier.iEntrance = iEntrance;
                stepIdentifier.stepId = stepId;
                readFollowingMap[orientedReadId].push_back(stepIdentifier);
            }
        }
    }

    for(uint64_t iExit=0; iExit<exits.size(); iExit++) {
        const auto& exit = exits[iExit];
        const AssemblyGraph::edge_descriptor e = exit.e;
        const AssemblyGraphEdge& edge = assemblyGraph[e];
        for(uint64_t stepId=0; stepId<edge.size(); stepId++) {
            const AssemblyGraphEdgeStep& step = edge[stepId];
            const AnchorPair& anchorPair = step.anchorPair;
            for(const OrientedReadId orientedReadId: anchorPair.orientedReadIds) {
                StepIdentifier stepIdentifier;
                stepIdentifier.iExit = iExit;
                stepIdentifier.stepId = stepId;
                readFollowingMap[orientedReadId].push_back(stepIdentifier);
            }
        }
    }

    if(debug) {
        for(const auto& p: readFollowingMap) {
            const OrientedReadId orientedReadId = p.first;
            const vector<StepIdentifier>& stepIdentifiers = p.second;
            cout << "Occurrences of " << orientedReadId << ":" << endl;
            for(const StepIdentifier& stepIdentifier: stepIdentifiers) {
                if(stepIdentifier.iExit == invalid<uint64_t>) {
                    cout << "In-" << stepIdentifier.iEntrance << "-" << stepIdentifier.stepId << " ";
                } else {
                    cout << "Out-" << stepIdentifier.iExit << "-" << stepIdentifier.stepId << " ";
                }
            }
            cout << endl;
        }
    }
}
