#include "TangleMatrix1.hpp"
#include "deduplicate.hpp"
using namespace shasta;

#include <iomanip>



TangleMatrix1::TangleMatrix1(
    const AssemblyGraph& assemblyGraph,
    vector<edge_descriptor> entrances,
    vector<edge_descriptor> exits,
    ostream& html) :
    assemblyGraph(assemblyGraph),
    entrances(entrances),
    exits(exits)
{
    // Sanity checks.
    SHASTA_ASSERT(not assemblyGraph.orientedReadSegments.empty());
    SHASTA_ASSERT(std::ranges::is_sorted(entrances, assemblyGraph.orderById));
    SHASTA_ASSERT(std::ranges::is_sorted(exits, assemblyGraph.orderById));

    // Gather the oriented reads that contribute to this TangleMatrix1.
    // These are the oriented reads that appear in at least one entrance
    // and one exit.
    gatherOrientedReads();
    if(html) {
        writeOrientedReads(html);
    }

    computeTotalTangleMatrix();
    if(html) {
        writeTotalTangleMatrix(html);
    }
}



void TangleMatrix1::gatherOrientedReads()
{
    // Find transitioning orientedReads that appear in one or more entrances.
    vector<OrientedReadId> entranceOrientedReadIds;
    for(const edge_descriptor e: entrances) {
        for(const auto& p: assemblyGraph[e].transitioningOrientedReadIds) {
            entranceOrientedReadIds.push_back(p.first);
        }
    }
    deduplicate(entranceOrientedReadIds);

    // Find transitioning orientedReads that appear in one or more exits.
    vector<OrientedReadId> exitOrientedReadIds;
    for(const edge_descriptor e: exits) {
        for(const auto& p: assemblyGraph[e].transitioningOrientedReadIds) {
            exitOrientedReadIds.push_back(p.first);
        }
    }
    deduplicate(exitOrientedReadIds);


    // Find OrientedReadIds that appear in at least one entrance
    // and at least one exit.
    vector<OrientedReadId> orientedReadIds;
    std::set_intersection(
        entranceOrientedReadIds.begin(), entranceOrientedReadIds.end(),
        exitOrientedReadIds.begin(), exitOrientedReadIds.end(),
        back_inserter(orientedReadIds));

    // Initialize the orientedReadInfos vector.
    for(const OrientedReadId orientedReadId: orientedReadIds) {
        orientedReadInfos.emplace_back(orientedReadId, entrances.size(), exits.size());
    }


    // Gather the number of appearances in entrances and exits for each of these oriented reads.
    for(uint64_t i=0; i<entrances.size(); i++) {
        const edge_descriptor e = entrances[i];
        for(const auto& p: assemblyGraph[e].transitioningOrientedReadIds) {
            const OrientedReadId orientedReadId = p.first;
            const uint64_t j = getOrientedReadIdIndex(orientedReadId);
            if(j != invalid<uint64_t>) {
                const uint64_t count = p.second;
                OrientedReadInfo& info = orientedReadInfos[j];
                info.entranceStepCount[i] = count;
            }
        }
    }
    for(uint64_t i=0; i<exits.size(); i++) {
        const edge_descriptor e = exits[i];
        for(const auto& p: assemblyGraph[e].transitioningOrientedReadIds) {
            const OrientedReadId orientedReadId = p.first;
            const uint64_t j = getOrientedReadIdIndex(orientedReadId);
            if(j != invalid<uint64_t>) {
                const uint64_t count = p.second;
                OrientedReadInfo& info = orientedReadInfos[j];
                info.exitStepCount[i] = count;
            }
        }
    }


    // For each of the oriented reads, compute the contribution to the
    // tangle matrix.
    for(OrientedReadInfo& orientedReadInfo: orientedReadInfos) {
        orientedReadInfo.computeTangleMatrix();
    }
}



TangleMatrix1::OrientedReadInfo::OrientedReadInfo(
    OrientedReadId orientedReadId,
    uint64_t entranceCount,
    uint64_t exitCount) :
    orientedReadId(orientedReadId),
    entranceStepCount(entranceCount, 0),
    exitStepCount(exitCount, 0)
{}



bool TangleMatrix1::OrientedReadInfo::operator<(const OrientedReadInfo& that) const
{
    return orientedReadId < that.orientedReadId;
}



void TangleMatrix1::OrientedReadInfo::computeTangleMatrix()
{
    const double entranceSum = double(std::accumulate(entranceStepCount.begin(), entranceStepCount.end(), 0UL));
    const double exitSum = double(std::accumulate(exitStepCount.begin(), exitStepCount.end(), 0UL));

    const uint64_t entranceCount = entranceStepCount.size();
    const uint64_t exitCount = exitStepCount.size();

    tangleMatrix.clear();
    tangleMatrix.resize(entranceCount, vector<double>(exitCount, 0));

    for(uint64_t i=0; i<entranceCount; i++) {
        for(uint64_t j=0; j<exitCount; j++) {
            tangleMatrix[i][j] = double(entranceStepCount[i] * exitStepCount[j]) / (entranceSum * exitSum);
        }
    }

}



// Return the index of a given OrientedReadId in the orientedReadInfos vector,
// or invalid<uint64_t> if not present.
uint64_t TangleMatrix1::getOrientedReadIdIndex(OrientedReadId orientedReadId) const
{
    const OrientedReadInfo target(orientedReadId);
    const auto it = std::lower_bound(orientedReadInfos.begin(), orientedReadInfos.end(), target);

    if((it == orientedReadInfos.end()) or (it->orientedReadId != orientedReadId)) {
        return invalid<uint64_t>;
    }

    return it - orientedReadInfos.begin();
}



void TangleMatrix1::writeOrientedReads(ostream& html) const
{
    html << std::setprecision(2) << std::defaultfloat << "<h3>Tangle matrix</h3>"
        "<h5>Oriented read contributions</h5>"
        "<table>"
        "<tr><th>Oriented<br>read<br>id";
    for(uint64_t i=0; i<entrances.size(); i++) {
        html << "<th>Entrance<br>" << assemblyGraph[entrances[i]].id;
    }
    for(uint64_t i=0; i<exits.size(); i++) {
        html << "<th>Exit<br>" << assemblyGraph[exits[i]].id;
    }
    html << "<th>Tangle<br>matrix";

    // Loop over oriented reads.
    for(const OrientedReadInfo& orientedReadInfo: orientedReadInfos) {
        const OrientedReadId orientedReadId = orientedReadInfo.orientedReadId;

        html << "<tr><th>" << orientedReadId;
        for(uint64_t i=0; i<entrances.size(); i++) {
            html << "<td class=centered>" << orientedReadInfo.entranceStepCount[i];
        }
        for(uint64_t j=0; j<exits.size(); j++) {
            html << "<td class=centered>" << orientedReadInfo.exitStepCount[j];
        }

        // Write the contribution of this oriented read to the total tangle matrix.
        html << "<td class=centered><table style='margin: 0 auto;'>";
        for(uint64_t i=0; i<entrances.size(); i++) {
            html << "<tr>";
            for(uint64_t j=0; j<exits.size(); j++) {
                html << "<td class=centered>" << orientedReadInfo.tangleMatrix[i][j];
            }
        }
        html << "</table>";
    }
    html << "</table>";
}



void TangleMatrix1::computeTotalTangleMatrix()
{
    const uint64_t entranceCount = entrances.size();
    const uint64_t exitCount = exits.size();

    tangleMatrix.clear();
    tangleMatrix.resize(entranceCount, vector<double>(exitCount, 0.));
    for(const OrientedReadInfo& orientedReadInfo: orientedReadInfos) {
        for(uint64_t i=0; i<entranceCount; i++) {
            for(uint64_t j=0; j<exitCount; j++) {
                tangleMatrix[i][j] += orientedReadInfo.tangleMatrix[i][j];
            }
        }
    }
}



void TangleMatrix1::writeTotalTangleMatrix(ostream& html) const
{
    html << "<h5>Total tangle matrix</h5><table><tr><th>";
    for(uint64_t j=0; j<exits.size(); j++) {
        html << "<th>" << assemblyGraph[exits[j]].id;
    }

    for(uint64_t i=0; i<entrances.size(); i++) {
        html << "<tr><th>" << assemblyGraph[entrances[i]].id;
        for(uint64_t j=0; j<exits.size(); j++) {
            html << "<td class=centered>" << tangleMatrix[i][j];
        }
    }

    html << "</table>";

}
