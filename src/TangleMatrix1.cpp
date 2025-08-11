#include "TangleMatrix1.hpp"
#include "deduplicate.hpp"
using namespace shasta;

#include <fstream.hpp>
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
    SHASTA_ASSERT(not assemblyGraph.orientedReadEdgeInformation.empty());
    SHASTA_ASSERT(std::ranges::is_sorted(entrances, assemblyGraph.orderById));
    SHASTA_ASSERT(std::ranges::is_sorted(exits, assemblyGraph.orderById));

    // Gather the oriented reads that contribute to this TangleMatrix1.
    // These are the oriented reads that appear in at least one entrance
    // and one exit.
    gatherCommonOrientedReads();
    if(html) {
        writeCommonOrientedReads(html);
    }

    computeTotalTangleMatrix();
    if(html) {
        writeTotalTangleMatrix(html);
    }
}



void TangleMatrix1::gatherCommonOrientedReads()
{
    const uint64_t entranceCount = entrances.size();
    const uint64_t exitCount = exits.size();

    // Gather AssemblyGraph::OrientedReadEdgeInfos for all entrances and exits.
    vector< vector<AssemblyGraph::OrientedReadEdgeInfo> > entranceInfos(entranceCount);
    for(uint64_t i=0; i<entranceCount; i++) {
        assemblyGraph.gatherOrientedReadInformationOnEdge(
            entrances[i],
            entranceInfos[i]);
    }
    vector< vector<AssemblyGraph::OrientedReadEdgeInfo> > exitInfos(exitCount);
    for(uint64_t i=0; i<exitCount; i++) {
        assemblyGraph.gatherOrientedReadInformationOnEdge(
            exits[i],
            exitInfos[i]);
    }

    // Find orientedReads that appear in one or more entrances.
    vector<OrientedReadId> entranceOrientedReadIds;
    for(const auto& v: entranceInfos) {
        for(const AssemblyGraph::OrientedReadEdgeInfo& info: v) {
            entranceOrientedReadIds.push_back(info.orientedReadId);
        }
    }
    deduplicate(entranceOrientedReadIds);

    // Find orientedReads that appear in one or more exits.
    vector<OrientedReadId> exitOrientedReadIds;
    for(const auto& v: exitInfos) {
        for(const AssemblyGraph::OrientedReadEdgeInfo& info: v) {
            exitOrientedReadIds.push_back(info.orientedReadId);
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



    // Loop over these common OrientedReadIds.
    // Maintain iterators into the entranceInfos and exitInfos.
    commonOrientedReadInfos.clear();
    vector< vector<AssemblyGraph::OrientedReadEdgeInfo>::const_iterator > entranceInfoIterators(entranceCount);
    for(uint64_t i=0; i<entranceCount; i++) {
        entranceInfoIterators[i] = entranceInfos[i].begin();
    }
    vector< vector<AssemblyGraph::OrientedReadEdgeInfo>::const_iterator > exitInfoIterators(exitCount);
    for(uint64_t i=0; i<exitCount; i++) {
        exitInfoIterators[i] = exitInfos[i].begin();
    }
    for(const OrientedReadId orientedReadId: orientedReadIds) {
        commonOrientedReadInfos.emplace_back(orientedReadId, entranceCount, exitCount);
        CommonOrientedReadInfo& commonOrientedReadInfo = commonOrientedReadInfos.back();

        for(uint64_t i=0; i<entranceCount; i++) {
            auto& it = entranceInfoIterators[i];
            while((it != entranceInfos[i].end()) and (it->orientedReadId < orientedReadId)) {
                ++it;
            }
            if((it == entranceInfos[i].end()) or (it->orientedReadId != orientedReadId)) {
                commonOrientedReadInfo.entranceStepCount[i] = 0;
            } else {
                commonOrientedReadInfo.entranceStepCount[i] = it->stepCount;
            }
        }

        for(uint64_t i=0; i<exitCount; i++) {
            auto& it = exitInfoIterators[i];
            while((it != exitInfos[i].end()) and (it->orientedReadId < orientedReadId)) {
                ++it;
            }
            if((it == exitInfos[i].end()) or (it->orientedReadId != orientedReadId)) {
                commonOrientedReadInfo.exitStepCount[i] = 0;
            } else {
                commonOrientedReadInfo.exitStepCount[i] = it->stepCount;
            }
        }

    }

    // For each of the oriented reads, compute the contribution to the
    // tangle matrix.
    for(CommonOrientedReadInfo& commonOrientedReadInfo: commonOrientedReadInfos) {
        commonOrientedReadInfo.computeTangleMatrix();
    }
}



TangleMatrix1::CommonOrientedReadInfo::CommonOrientedReadInfo(
    OrientedReadId orientedReadId,
    uint64_t entranceCount,
    uint64_t exitCount) :
    orientedReadId(orientedReadId),
    entranceStepCount(entranceCount, 0),
    exitStepCount(exitCount, 0)
{}



bool TangleMatrix1::CommonOrientedReadInfo::operator<(const CommonOrientedReadInfo& that) const
{
    return orientedReadId < that.orientedReadId;
}



void TangleMatrix1::CommonOrientedReadInfo::computeTangleMatrix()
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



// Return the index of a given OrientedReadId in the commonOrientedReadInfos vector,
// or invalid<uint64_t> if not present.
uint64_t TangleMatrix1::getCommonOrientedReadIdIndex(OrientedReadId orientedReadId) const
{
    const CommonOrientedReadInfo target(orientedReadId);
    const auto it = std::lower_bound(commonOrientedReadInfos.begin(), commonOrientedReadInfos.end(), target);

    if((it == commonOrientedReadInfos.end()) or (it->orientedReadId != orientedReadId)) {
        return invalid<uint64_t>;
    }

    return it - commonOrientedReadInfos.begin();
}



void TangleMatrix1::writeCommonOrientedReads(ostream& html) const
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

    // Loop over common oriented reads.
    for(const CommonOrientedReadInfo& commonOrientedReadInfo: commonOrientedReadInfos) {
        const OrientedReadId orientedReadId = commonOrientedReadInfo.orientedReadId;

        html << "<tr><th>" << orientedReadId;
        for(uint64_t i=0; i<entrances.size(); i++) {
            html << "<td class=centered>" << commonOrientedReadInfo.entranceStepCount[i];
        }
        for(uint64_t j=0; j<exits.size(); j++) {
            html << "<td class=centered>" << commonOrientedReadInfo.exitStepCount[j];
        }

        // Write the contribution of this oriented read to the total tangle matrix.
        html << "<td class=centered><table style='margin: 0 auto;'>";
        for(uint64_t i=0; i<entrances.size(); i++) {
            html << "<tr>";
            for(uint64_t j=0; j<exits.size(); j++) {
                html << "<td class=centered>" << commonOrientedReadInfo.tangleMatrix[i][j];
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
    for(const CommonOrientedReadInfo& commonOrientedReadInfo: commonOrientedReadInfos) {
        for(uint64_t i=0; i<entranceCount; i++) {
            for(uint64_t j=0; j<exitCount; j++) {
                tangleMatrix[i][j] += commonOrientedReadInfo.tangleMatrix[i][j];
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
