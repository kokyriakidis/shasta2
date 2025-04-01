#include "TangleMatrix.hpp"
#include "Anchor.hpp"
#include "deduplicate.hpp"
using namespace shasta;



TangleMatrix::TangleMatrix(
    const AssemblyGraph& assemblyGraph,
    vector<edge_descriptor> entranceEdges,
    vector<edge_descriptor> exitEdges)
{
    construct(assemblyGraph, entranceEdges, exitEdges);
}



void TangleMatrix::construct(
    const AssemblyGraph& assemblyGraph,
    vector<edge_descriptor> entranceEdges,
    vector<edge_descriptor> exitEdges)
{
    // Sort the edge descriptors by id.
    AssemblyGraphEdgeOrderById sorter(assemblyGraph);
    sort(entranceEdges.begin(), entranceEdges.end(), sorter);
    sort(exitEdges.begin(), exitEdges.end(), sorter);

    // Check for duplicates.
    AssemblyGraphEdgeCompareEqualById comparator(assemblyGraph);
    SHASTA_ASSERT(std::adjacent_find(
        entranceEdges.begin(), entranceEdges.end(), comparator) == entranceEdges.end());
    SHASTA_ASSERT(std::adjacent_find(
        exitEdges.begin(), exitEdges.end(), comparator) == exitEdges.end());

    // Create the entrances.
    // Store the second to last AnchorId of each entrance.
    for(const edge_descriptor e: entranceEdges) {
        const AssemblyGraphEdge& edge = assemblyGraph[e];
        const AnchorId anchorId = edge.secondToLast();
        const uint64_t coverage = assemblyGraph.anchors[anchorId].coverage();
        entrances.emplace_back(e, anchorId, coverage);
    }

    // Create the exits.
    // Store the second AnchorId of each exit.
    for(const edge_descriptor e: exitEdges) {
        const AssemblyGraphEdge& edge = assemblyGraph[e];
        SHASTA_ASSERT(not edge.empty());
        const AnchorId anchorId = edge.second();
        const uint64_t coverage = assemblyGraph.anchors[anchorId].coverage();
        exits.emplace_back(e, anchorId, coverage);
    }


    // Now we can compute the tangle matrix and store common coverage
    // for each entrance and exit.
    tangleMatrix.resize(entrances.size(), vector<uint64_t>(exits.size(), 0));
    for(uint64_t iEntrance=0; iEntrance<entrances.size(); iEntrance++) {
        Entrance& entrance = entrances[iEntrance];
        for(uint64_t iExit=0; iExit<exits.size(); iExit++) {
            Exit& exit = exits[iExit];
            const uint64_t commonCount = assemblyGraph.anchors.countCommon(entrance.anchorId, exit.anchorId);
            tangleMatrix[iEntrance][iExit] = commonCount;
            entrance.commonCoverage += commonCount;
            exit.commonCoverage += commonCount;
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
        "<tr><td rowspan=6 colspan=6>"
        "<th colspan=" << (exits.size() + 1) << " style='background-color:LightCyan'>Exits";

    html << "<tr><th style='background-color:LightCyan'>Index";
    for(uint64_t i=0; i<exits.size(); i++) {
        html << "<td class=centered style='background-color:LightCyan'>" << i;
    }

    html << "<tr><th style='background-color:LightCyan'>Segment";
    for(uint64_t i=0; i<exits.size(); i++) {
        html << "<td class=centered style='background-color:LightCyan'>" << assemblyGraph[exits[i].e].id;
    }

    html << "<tr><th style='background-color:LightCyan'>AnchorId";
    for(uint64_t i=0; i<exits.size(); i++) {
        html << "<td class=centered style='background-color:LightCyan'>" << anchorIdToString(exits[i].anchorId);
    }

    html << "<tr><th style='background-color:LightCyan'>Coverage";
    for(uint64_t i=0; i<exits.size(); i++) {
        html << "<td class=centered style='background-color:LightCyan'>" << exits[i].coverage;
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
        "<th style='background-color:CornSilk'>AnchorId"
        "<th style='background-color:CornSilk'>Coverage"
        "<th style='background-color:CornSilk'>Common<br>coverage<td>";

    for(uint64_t i=0; i<entrances.size(); i++) {
        html << "<td>";
    }

    for(uint64_t iEntrance=0; iEntrance<entrances.size(); iEntrance++) {
        const Entrance& entrance = entrances[iEntrance];
        html <<
            "<tr>"
            "<td class=centered style='background-color:CornSilk'>" << iEntrance <<
            "<td class=centered style='background-color:CornSilk'>" << assemblyGraph[entrance.e].id <<
            "<td class=centered style='background-color:CornSilk'>" << anchorIdToString(entrance.anchorId) <<
            "<td class=centered style='background-color:CornSilk'>" << entrance.coverage <<
            "<td class=centered style='background-color:CornSilk'>" << entrance.commonCoverage <<
            "<td>";
        for(uint64_t iExit=0; iExit<exits.size(); iExit++) {
            html << "<td class=centered style='background-color:LightPink'>" << tangleMatrix[iEntrance][iExit];
        }
    }

    html << "</table>";
}

