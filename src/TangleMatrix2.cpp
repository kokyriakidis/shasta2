#include "TangleMatrix2.hpp"
#include "Anchor.hpp"
#include "Base.hpp"
using namespace shasta;



TangleMatrix2::TangleMatrix2(
    const AssemblyGraph2& assemblyGraph2,
    vector<vertex_descriptor> entranceVertices,
    vector<vertex_descriptor> exitVertices,
    double aDrift,
    double bDrift)
{

    for(const vertex_descriptor v: entranceVertices) {
        const AssemblyGraph2Vertex& vertex = assemblyGraph2[v];
        entrances.emplace_back(v, vertex.back());
    }

    for(const vertex_descriptor v: exitVertices) {
        const AssemblyGraph2Vertex& vertex = assemblyGraph2[v];
        exits.emplace_back(v, vertex.front());
    }

    // Check that no entrance/exit pairs are adjacent.
    for(const auto& entrance: entrances) {
        for(const auto& exit: exits) {
            if(entrance.step.anchorPair.anchorIdB == exit.step.anchorPair.anchorIdA) {
                throw runtime_error(
                    "Entrance " + to_string(assemblyGraph2[entrance.v].id) +
                    " is adjacent to exit " + to_string(assemblyGraph2[exit.v].id));
            }
        }
    }


    // Compute the tangle matrix.
    tangleMatrix.resize(entrances.size(), vector<AnchorPair>(exits.size()));
    for(uint64_t iEntrance=0; iEntrance<entrances.size(); iEntrance++) {
        EntranceOrExit& entrance = entrances[iEntrance];
        for(uint64_t iExit=0; iExit<exits.size(); iExit++) {
            EntranceOrExit& exit = exits[iExit];

            tangleMatrix[iEntrance][iExit] = assemblyGraph2.anchors.bridge(
                entrance.step.anchorPair,
                exit.step.anchorPair,
                aDrift, bDrift);
            const uint64_t coverage = tangleMatrix[iEntrance][iExit].orientedReadIds.size();
            entrance.commonCoverage += coverage;
            exit.commonCoverage += coverage;
        }
    }
}



void TangleMatrix2::writeHtml(
    const AssemblyGraph2& assemblyGraph2,
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
        html << "<td class=centered style='background-color:LightCyan'>" << assemblyGraph2[exits[i].v].id;
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
            "<td class=centered style='background-color:CornSilk'>" << assemblyGraph2[entrance.v].id <<
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
