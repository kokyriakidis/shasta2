#include "TangleMatrix3.hpp"
#include "Anchor.hpp"
#include "Base.hpp"
using namespace shasta;



TangleMatrix3::TangleMatrix3(
    const AssemblyGraph& assemblyGraph3,
    vector<edge_descriptor> entranceEdges,
    vector<edge_descriptor> exitEdges,
    double aDrift,
    double bDrift)
{

    for(const edge_descriptor e: entranceEdges) {
        const AssemblyGraphEdge& edge = assemblyGraph3[e];
        entrances.emplace_back(e, edge.back());
    }

    for(const edge_descriptor e: exitEdges) {
        const AssemblyGraphEdge& edge = assemblyGraph3[e];
        exits.emplace_back(e, edge.front());
    }



    // Compute the tangle matrix.
    tangleMatrix.resize(entrances.size(), vector<AnchorPair>(exits.size()));
    for(uint64_t iEntrance=0; iEntrance<entrances.size(); iEntrance++) {
        EntranceOrExit& entrance = entrances[iEntrance];
        for(uint64_t iExit=0; iExit<exits.size(); iExit++) {
            EntranceOrExit& exit = exits[iExit];

            tangleMatrix[iEntrance][iExit] = assemblyGraph3.anchors.bridge(
                entrance.step.anchorPair,
                exit.step.anchorPair,
                aDrift, bDrift);
            const uint64_t coverage = tangleMatrix[iEntrance][iExit].orientedReadIds.size();
            entrance.commonCoverage += coverage;
            exit.commonCoverage += coverage;
        }
    }
}



void TangleMatrix3::writeHtml(
    const AssemblyGraph& assemblyGraph3,
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
        html << "<td class=centered style='background-color:LightCyan'>" << assemblyGraph3[exits[i].e].id;
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
            "<td class=centered style='background-color:CornSilk'>" << assemblyGraph3[entrance.e].id <<
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

