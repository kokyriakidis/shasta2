#if 0
#include "TangleMatrix.hpp"
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

    // Find OrientedReadIds that appear in more than one entrance.
    // These will be discarded when creating the Entrances.
    vector<uint64_t> count;
    for(const edge_descriptor e: entranceEdges) {
        const AssemblyGraphEdge& edge = assemblyGraph[e];
        SHASTA_ASSERT(not edge.empty());
        const vector<OrientedReadId>& orientedReadIds = edge.back().orientedReadIds;
        copy(orientedReadIds.begin(), orientedReadIds.end(),
            back_inserter(duplicateOrientedReadIdsOnEntrances));
    }
    deduplicateAndCountWithThreshold(duplicateOrientedReadIdsOnEntrances, count, 2UL);

    // Find OrientedReadIds that appear in more than one exit.
    // These will be discarded when creating the Exits.
    for(const edge_descriptor e: exitEdges) {
        const AssemblyGraphEdge& edge = assemblyGraph[e];
        SHASTA_ASSERT(not edge.empty());
        const vector<OrientedReadId>& orientedReadIds = edge.front().orientedReadIds;
        copy(orientedReadIds.begin(), orientedReadIds.end(),
            back_inserter(duplicateOrientedReadIdsOnExits));
    }
    deduplicateAndCountWithThreshold(duplicateOrientedReadIdsOnExits, count, 2UL);



    // Create the entrances.
    // Store the last AssemblyGraphStep of each entrance.
    for(const edge_descriptor e: entranceEdges) {
        const AssemblyGraphEdge& edge = assemblyGraph[e];
        SHASTA_ASSERT(not edge.empty());
        entrances.emplace_back(e, edge.back(), duplicateOrientedReadIdsOnEntrances);
    }

    // Create the exits.
    // Store the first AssemblyGraphStep of each exit.
    for(const edge_descriptor e: exitEdges) {
        const AssemblyGraphEdge& edge = assemblyGraph[e];
        SHASTA_ASSERT(not edge.empty());
        exits.emplace_back(e, edge.front(), duplicateOrientedReadIdsOnExits);
    }

    // Create the AnchorPairs that we would get if we were to join
    // each Entrance with each exit.
    // Some of these will be used when detangling.
    joinedAnchorPairs.resize(entrances.size(), vector<AnchorPair>(exits.size()));
    for(uint64_t iEntrance=0; iEntrance<entrances.size(); iEntrance++) {
        const AnchorPair& entranceAnchorPair = entrances[iEntrance].step;
        for(uint64_t iExit=0; iExit<exits.size(); iExit++) {
            const AnchorPair& exitAnchorPair = exits[iExit].step;
            joinedAnchorPairs[iEntrance][iExit] = AnchorPair(assemblyGraph.anchors,
                entranceAnchorPair, exitAnchorPair);
        }
    }

    // Now we can compute the tangle matrix.
    vector<OrientedReadId> commonOrientedReadIds;
    tangleMatrix.resize(entrances.size(), vector<uint64_t>(exits.size(), 0));
    for(uint64_t iEntrance=0; iEntrance<entrances.size(); iEntrance++) {
        for(uint64_t iExit=0; iExit<exits.size(); iExit++) {
            tangleMatrix[iEntrance][iExit] = joinedAnchorPairs[iEntrance][iExit].size();
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
        "<tr><th>Index<th>Segment<th>Coverage<th>Usable<br>coverage";
    for(uint64_t iEntrance=0; iEntrance<entrances.size(); iEntrance++) {
        const Entrance& entrance = entrances[iEntrance];
        const edge_descriptor e = entrance.e;
        const AssemblyGraphEdge& edge = assemblyGraph[e];
        const AssemblyGraphStep& edgeLastStep = edge.back();
        html <<
            "<tr>"
            "<td class = centered>" << iEntrance <<
            "<td class = centered>" << assemblyGraph[entrance.e].id <<
            "<td class = centered>" << edgeLastStep.orientedReadIds.size() <<
            "<td class = centered>" << entrance.step.orientedReadIds.size();
    }
    html << "</table>";



    // Exits.
    html <<
        "<h3>Exits</h3>"
        "<p><table>"
        "<tr><th>Index<th>Segment<th>Coverage<th>Usable<br>coverage";
    for(uint64_t iExit=0; iExit<exits.size(); iExit++) {
        const Exit& exit = exits[iExit];
        const edge_descriptor e = exit.e;
        const AssemblyGraphEdge& edge = assemblyGraph[e];
        const AssemblyGraphStep& edgeFirstStep = edge.front();
        html <<
            "<tr>"
            "<td class = centered>" << iExit <<
            "<td class = centered>" << assemblyGraph[exit.e].id <<
            "<td class = centered>" << edgeFirstStep.orientedReadIds.size() <<
            "<td class = centered>" << exit.step.orientedReadIds.size();
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
    html << "<t/able>";
}
#endif
