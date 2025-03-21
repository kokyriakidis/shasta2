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

    // Create the entrances.
    // Store the last AssemblyGraphStep of each entrance.
    for(const edge_descriptor e: entranceEdges) {
        const AssemblyGraphEdge& edge = assemblyGraph[e];
        SHASTA_ASSERT(not edge.empty());
        entrances.emplace_back(e, edge.back());
    }

    // Create the exits.
    // Store the first AssemblyGraphStep of each exit.
    for(const edge_descriptor e: exitEdges) {
        const AssemblyGraphEdge& edge = assemblyGraph[e];
        SHASTA_ASSERT(not edge.empty());
        exits.emplace_back(e, edge.front());
    }

    // If an OrientedReadId appears in more than one Entrance,
    // we ignore it when creating the TangleMatrix.
    vector<uint64_t> count;
    for(const Entrance& entrance: entrances) {
        const vector<OrientedReadId>& orientedReadIds = entrance.step.orientedReadIds;
        copy(orientedReadIds.begin(), orientedReadIds.end(),
            back_inserter(duplicateOrientedReadIdsOnEntrances));
    }
    deduplicateAndCountWithThreshold(duplicateOrientedReadIdsOnEntrances, count, 2UL);
    duplicateOrientedReadIdsOnEntrances.shrink_to_fit();

    // If an OrientedReadId appears in more than one Exit,
    // we ignore it when creating the TangleMatrix.
    for(const Exit& exit: exits) {
        const vector<OrientedReadId>& orientedReadIds = exit.step.orientedReadIds;
        copy(orientedReadIds.begin(), orientedReadIds.end(),
            back_inserter(duplicateOrientedReadIdsOnExits));
    }
    deduplicateAndCountWithThreshold(duplicateOrientedReadIdsOnExits, count, 2UL);
    duplicateOrientedReadIdsOnExits.shrink_to_fit();

    // Now we can fill in the OrientedReadIds of each Entrance and Exit.
    for(Entrance& entrance: entrances) {
        std::set_difference(
            entrance.step.orientedReadIds.begin(), entrance.step.orientedReadIds.end(),
            duplicateOrientedReadIdsOnEntrances.begin(), duplicateOrientedReadIdsOnEntrances.end(),
            back_inserter(entrance.orientedReadIds)
            );
    }
    for(Exit& exit: exits) {
        std::set_difference(
            exit.step.orientedReadIds.begin(), exit.step.orientedReadIds.end(),
            duplicateOrientedReadIdsOnExits.begin(), duplicateOrientedReadIdsOnExits.end(),
            back_inserter(exit.orientedReadIds)
            );
    }

    // Now we can compute the tangle matrix.
    vector<OrientedReadId> commonOrientedReadIds;
    tangleMatrix.resize(entrances.size(), vector<uint64_t>(exits.size(), 0));
    for(uint64_t iEntrance=0; iEntrance<entrances.size(); iEntrance++) {
        const vector<OrientedReadId>& entranceOrientedReadIds = entrances[iEntrance].orientedReadIds;
        for(uint64_t iExit=0; iExit<exits.size(); iExit++) {
            const vector<OrientedReadId>& exitOrientedReadIds = exits[iExit].orientedReadIds;
            commonOrientedReadIds.clear();
            std::set_intersection(
                entranceOrientedReadIds.begin(), entranceOrientedReadIds.end(),
                exitOrientedReadIds.begin(), exitOrientedReadIds.end(),
                back_inserter(commonOrientedReadIds));
            tangleMatrix[iEntrance][iExit] = commonOrientedReadIds.size();
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
        html <<
            "<tr>"
            "<td class = centered>" << iEntrance <<
            "<td class = centered>" << assemblyGraph[entrance.e].id <<
            "<td class = centered>" << entrance.step.size() <<
            "<td class = centered>" << entrance.orientedReadIds.size();
    }
    html << "</table>";



    // Exits.
    html <<
        "<h3>Exits</h3>"
        "<p><table>"
        "<tr><th>Index<th>Segment<th>Coverage<th>Usable<br>coverage";
    for(uint64_t iExit=0; iExit<exits.size(); iExit++) {
        const Exit& exit = exits[iExit];
        html <<
            "<tr>"
            "<td class = centered>" << iExit <<
            "<td class = centered>" << assemblyGraph[exit.e].id <<
            "<td class = centered>" << exit.step.size() <<
            "<td class = centered>" << exit.orientedReadIds.size();
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
