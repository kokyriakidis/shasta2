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
    html << "<h2>Tangle matrix</h2>";



    // Entrances.
    html <<
        "<h3>Entrances</h3>"
        "<p><table>"
        "<tr><th>Index<th>Segment<th>AnchorId<th>Coverage";
    for(uint64_t iEntrance=0; iEntrance<entrances.size(); iEntrance++) {
        const Entrance& entrance = entrances[iEntrance];
        const edge_descriptor e = entrance.e;
        const AnchorId anchorId = entrance.anchorId;
        const uint64_t coverage = assemblyGraph.anchors[anchorId].coverage();
        html <<
            "<tr>"
            "<td class = centered>" << iEntrance <<
            "<td class = centered>" << assemblyGraph[e].id <<
            "<td class = centered>" << anchorIdToString(anchorId) <<
            "<td class = centered>" << coverage;
    }
    html << "</table>";



    // Exits.
    html <<
        "<h3>Exits</h3>"
        "<p><table>"
        "<tr><th>Index<th>Segment<th>AnchorId<th>Coverage";
    for(uint64_t iExit=0; iExit<exits.size(); iExit++) {
        const Exit& exit = exits[iExit];
        const edge_descriptor e = exit.e;
        const AnchorId anchorId = exit.anchorId;
        const uint64_t coverage = assemblyGraph.anchors[anchorId].coverage();
        html <<
            "<tr>"
            "<td class = centered>" << iExit <<
            "<td class = centered>" << assemblyGraph[e].id <<
            "<td class = centered>" << anchorIdToString(anchorId) <<
            "<td class = centered>" << coverage;
    }
    html << "</table>";



    // Tangle matrix.
    html <<
        "<h3>Tangle matrix</h3>"
        "<p>One row for each entrance, one column for each exit."
        "<p><table>"
        "<tr><td rowspan=2 colspan=2>";
    for(uint64_t iExit=0; iExit<exits.size(); iExit++) {
        html << "<th class=centered>" << iExit;
    }
    html << "<tr>";
    for(uint64_t iExit=0; iExit<exits.size(); iExit++) {
        const Exit& exit = exits[iExit];
        html << "<th class=centered>" << assemblyGraph[exit.e].id;
    }
    for(uint64_t iEntrance=0; iEntrance<entrances.size(); iEntrance++) {
        const Entrance& entrance = entrances[iEntrance];
        html <<
            "<tr>"
            "<th class=centered>" << iEntrance <<
            "<th class=centered>" << assemblyGraph[entrance.e].id;
        for(uint64_t iExit=0; iExit<exits.size(); iExit++) {
            html << "<td class=centered>" << tangleMatrix[iEntrance][iExit];
        }
    }
    html << "</table>";
}

