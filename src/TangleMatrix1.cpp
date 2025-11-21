#include "TangleMatrix1.hpp"
#include "Anchor.hpp"
#include "deduplicate.hpp"
#include "Options.hpp"
using namespace shasta;

#include "fstream.hpp"
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
    SHASTA2_ASSERT(std::ranges::is_sorted(entrances, assemblyGraph.orderById));
    SHASTA2_ASSERT(std::ranges::is_sorted(exits, assemblyGraph.orderById));

    // Gather oriented reads in the representative region of each entrance and exit.
    gatherOrientedReads(assemblyGraph.options.representativeRegionStepCount);
    if(html) {
        html << std::setprecision(2) << std::defaultfloat << "<h3>Tangle matrix</h3>";
        writeOrientedReads(html);
    }

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



void TangleMatrix1::gatherOrientedReads(uint64_t representativeRegionLength)
{
    const uint64_t entranceCount = entrances.size();
    const uint64_t exitCount = exits.size();

    entranceOrientedReadInfos.resize(entranceCount);
    for(uint64_t iEntrance=0; iEntrance<entranceCount; iEntrance++) {
        gatherEntranceOrientedReads(iEntrance, representativeRegionLength);
    }
    exitOrientedReadInfos.resize(exitCount);
    for(uint64_t iExit=0; iExit<exitCount; iExit++) {
        gatherExitOrientedReads(iExit, representativeRegionLength);
    }

}



void TangleMatrix1::gatherEntranceOrientedReads(
    uint64_t iEntrance,
    uint64_t representativeRegionLength)
{
    // Get the edge for this entrance.
    const edge_descriptor e = entrances[iEntrance];
    const AssemblyGraphEdge& edge = assemblyGraph[e];

    // Define the beginning of the representative region.
    uint64_t begin = invalid<uint64_t>;
    const uint64_t stepCount = edge.size();
    if(representativeRegionLength >= stepCount) {
        begin = 0;
    } else {
        begin = stepCount - representativeRegionLength;
    }



    // A work vector to store information on each appearance of an OrientedRead
    // on a step of the representative region.
    class WorkInfo {
    public:
        OrientedReadId orientedReadId;
        uint32_t positionInJourney;
        WorkInfo(
            OrientedReadId orientedReadId,
            uint32_t positionInJourney) :
            orientedReadId(orientedReadId),
            positionInJourney(positionInJourney)
        {}
        bool operator<(const WorkInfo& that) const
        {
            return orientedReadId < that.orientedReadId;
        }
    };
    vector<WorkInfo> workInfos;


    // Loop over the representative region and fill in the workInfos vector.
    for(uint64_t stepId=begin; stepId<stepCount; stepId++) {
        const AssemblyGraphEdgeStep& step = edge[stepId];
        const AnchorPair& anchorPair = step.anchorPair;
        const AnchorId anchorIdB = anchorPair.anchorIdB;
        const Anchor anchorB = assemblyGraph.anchors[anchorIdB];

        // Joint loop over the OrientedReadIds in the AnchorPair
        // and the OrientedReadIds in anchorB.
        auto itB = anchorB.begin();
        for(const OrientedReadId orientedReadId: anchorPair.orientedReadIds) {
            while(itB->orientedReadId < orientedReadId) {
                ++itB;
                SHASTA2_ASSERT(itB != anchorB.end());
            }
            SHASTA2_ASSERT(itB->orientedReadId == orientedReadId);
            workInfos.emplace_back(orientedReadId, itB->positionInJourney);
        }
    }

    // Sort the WorkInfos by OrientedReadId.
    sort(workInfos.begin(), workInfos.end());



    // Each streak with the same OrientedReadId generates an OrientedReadInfo.
    vector<OrientedReadInfo>& orientedReadInfos = entranceOrientedReadInfos[iEntrance];
    for(auto streakBegin=workInfos.begin(); streakBegin!=workInfos.end(); /* Update later */) {
        const OrientedReadId orientedReadId = streakBegin->orientedReadId;

        auto streakEnd = streakBegin;
        while((streakEnd!= workInfos.end()) and (streakEnd->orientedReadId == orientedReadId)) {
            ++streakEnd;
        }

        // Loop over this streak.
        OrientedReadInfo orientedReadInfo;
        orientedReadInfo.orientedReadId = orientedReadId;
        orientedReadInfo.stepCount = streakEnd - streakBegin;
        orientedReadInfo.positionInJourney = 0;
        for(auto it=streakBegin; it!=streakEnd; it++) {
            orientedReadInfo.positionInJourney = max(orientedReadInfo.positionInJourney, it->positionInJourney);
        }
        orientedReadInfos.push_back(orientedReadInfo);

        // Prepare to process the next streak.
        streakBegin = streakEnd;
    }
}



void TangleMatrix1::gatherExitOrientedReads(
    uint64_t iExit,
    uint64_t representativeRegionLength)
{
    // Get the edge for this entrance.
    const edge_descriptor e = exits[iExit];
    const AssemblyGraphEdge& edge = assemblyGraph[e];

    // Define the end of the representative region.
    const uint64_t stepCount = edge.size();
    const uint64_t end = min(stepCount, representativeRegionLength);




    // A work vector to store information on each appearance of an OrientedRead
    // on a step of the representative region.
    class WorkInfo {
    public:
        OrientedReadId orientedReadId;
        uint32_t positionInJourney;
        WorkInfo(
            OrientedReadId orientedReadId,
            uint32_t positionInJourney) :
            orientedReadId(orientedReadId),
            positionInJourney(positionInJourney)
        {}
        bool operator<(const WorkInfo& that) const
        {
            return orientedReadId < that.orientedReadId;
        }
    };
    vector<WorkInfo> workInfos;


    // Loop over the representative region and fill in the workInfos vector.
    for(uint64_t stepId=0; stepId<end; stepId++) {
        const AssemblyGraphEdgeStep& step = edge[stepId];
        const AnchorPair& anchorPair = step.anchorPair;
        const AnchorId anchorIdA = anchorPair.anchorIdA;
        const Anchor anchorA = assemblyGraph.anchors[anchorIdA];

        // Joint loop over the OrientedReadIds in the AnchorPair
        // and the OrientedReadIds in anchorA.
        auto itA = anchorA.begin();
        for(const OrientedReadId orientedReadId: anchorPair.orientedReadIds) {
            while(itA->orientedReadId < orientedReadId) {
                ++itA;
                SHASTA2_ASSERT(itA != anchorA.end());
            }
            SHASTA2_ASSERT(itA->orientedReadId == orientedReadId);
            workInfos.emplace_back(orientedReadId, itA->positionInJourney);
        }
    }

    // Sort the WorkInfos by OrientedReadId.
    sort(workInfos.begin(), workInfos.end());



    // Each streak with the same OrientedReadId generates an OrientedReadInfo.
    vector<OrientedReadInfo>& orientedReadInfos = exitOrientedReadInfos[iExit];
    for(auto streakBegin=workInfos.begin(); streakBegin!=workInfos.end(); /* Update later */) {
        const OrientedReadId orientedReadId = streakBegin->orientedReadId;

        auto streakEnd = streakBegin;
        while((streakEnd!= workInfos.end()) and (streakEnd->orientedReadId == orientedReadId)) {
            ++streakEnd;
        }

        // Loop over this streak.
        OrientedReadInfo orientedReadInfo;
        orientedReadInfo.orientedReadId = orientedReadId;
        orientedReadInfo.stepCount = streakEnd - streakBegin;
        orientedReadInfo.positionInJourney = std::numeric_limits<uint32_t>::max();
        for(auto it=streakBegin; it!=streakEnd; it++) {
            orientedReadInfo.positionInJourney = min(orientedReadInfo.positionInJourney, it->positionInJourney);
        }
        orientedReadInfos.push_back(orientedReadInfo);

        // Prepare to process the next streak.
        streakBegin = streakEnd;
    }
}


void TangleMatrix1::writeOrientedReads(ostream& html) const
{
    const uint64_t entranceCount = entrances.size();
    const uint64_t exitCount = exits.size();

    for(uint64_t iEntrance=0; iEntrance<entranceCount; iEntrance++) {
        const edge_descriptor e = entrances[iEntrance];
        html << "<h4>Oriented reads on entrance " << assemblyGraph[e].id << "</h4>"
            "These are the oriented reads that appear in the " <<
            "representative region of this entrance."
            "<p><table><tr>"
            "<th>Oriented<br>read<br>Id"
            "<th>Number<br>of<br>steps"
            "<th>Max<br>position<br>in journey";
        for(const OrientedReadInfo& orientedReadInfo: entranceOrientedReadInfos[iEntrance]) {
            html <<
                "<tr>"
                "<td class=centered>" << orientedReadInfo.orientedReadId <<
                "<td class=centered>" << orientedReadInfo.stepCount <<
                "<td class=centered>" << orientedReadInfo.positionInJourney;
        }
        html << "</table>";
    }

    for(uint64_t iExit=0; iExit<exitCount; iExit++) {
        const edge_descriptor e = exits[iExit];
        html << "<h4>Oriented reads on exit " << assemblyGraph[e].id << "</h4>"
            "<p>These are the oriented reads that appear in the " <<
            "representative region of this exit."
            "<table><tr>"
            "<th>Oriented<br>read<br>Id"
            "<th>Number<br>of<br>steps"
            "<th>Min<br>position<br>in journey";
        for(const OrientedReadInfo& orientedReadInfo: exitOrientedReadInfos[iExit]) {
            html <<
                "<tr>"
                "<td class=centered>" << orientedReadInfo.orientedReadId <<
                "<td class=centered>" << orientedReadInfo.stepCount <<
                "<td class=centered>" << orientedReadInfo.positionInJourney;
        }
        html << "</table>";
    }

}


#if 0
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
#endif



void TangleMatrix1::gatherCommonOrientedReads()
{
    const uint64_t entranceCount = entrances.size();
    const uint64_t exitCount = exits.size();

    // Find orientedReads that appear in one or more entrances.
    vector<OrientedReadId> entranceOrientedReadIds;
    for(const vector<OrientedReadInfo>& orientedReadInfos: entranceOrientedReadInfos) {
        for(const OrientedReadInfo& orientedReadInfo: orientedReadInfos) {
            entranceOrientedReadIds.push_back(orientedReadInfo.orientedReadId);
        }
    }
    deduplicate(entranceOrientedReadIds);

    // Find orientedReads that appear in one or more entrances.
    vector<OrientedReadId> exitOrientedReadIds;
    for(const vector<OrientedReadInfo>& orientedReadInfos: exitOrientedReadInfos) {
        for(const OrientedReadInfo& orientedReadInfo: orientedReadInfos) {
            exitOrientedReadIds.push_back(orientedReadInfo.orientedReadId);
        }
    }
    deduplicate(exitOrientedReadIds);

    // Find OrientedReadIds that appear in at least one entrance
    // and at least one exit.
    vector<OrientedReadId> commonOrientedReadIds;
    std::set_intersection(
        entranceOrientedReadIds.begin(), entranceOrientedReadIds.end(),
        exitOrientedReadIds.begin(), exitOrientedReadIds.end(),
        back_inserter(commonOrientedReadIds));



    // Loop over these common OrientedReadIds.
    // Maintain iterators into the entranceOrientedReadInfos and exitOrientedReadInfos.
    commonOrientedReadInfos.clear();
    vector< vector<OrientedReadInfo>::const_iterator > entranceInfoIterators(entranceCount);
    for(uint64_t i=0; i<entranceCount; i++) {
        entranceInfoIterators[i] = entranceOrientedReadInfos[i].begin();
    }
    vector< vector<OrientedReadInfo>::const_iterator > exitInfoIterators(exitCount);
    for(uint64_t i=0; i<exitCount; i++) {
        exitInfoIterators[i] = exitOrientedReadInfos[i].begin();
    }

    for(const OrientedReadId orientedReadId: commonOrientedReadIds) {
        commonOrientedReadInfos.emplace_back(orientedReadId, entranceCount, exitCount);
        CommonOrientedReadInfo& commonOrientedReadInfo = commonOrientedReadInfos.back();
        commonOrientedReadInfo.maxPositionInJourneyOnEntrances = 0;
        commonOrientedReadInfo.minPositionInJourneyOnExits = std::numeric_limits<uint32_t>::max();

        for(uint64_t i=0; i<entranceCount; i++) {
            auto& it = entranceInfoIterators[i];
            while((it != entranceOrientedReadInfos[i].end()) and (it->orientedReadId < orientedReadId)) {
                ++it;
            }
            if((it == entranceOrientedReadInfos[i].end()) or (it->orientedReadId != orientedReadId)) {
                commonOrientedReadInfo.entranceStepCount[i] = 0;
            } else {
                commonOrientedReadInfo.entranceStepCount[i] = it->stepCount;
                commonOrientedReadInfo.maxPositionInJourneyOnEntrances =
                    max(commonOrientedReadInfo.maxPositionInJourneyOnEntrances, it->positionInJourney);
            }
        }

        for(uint64_t i=0; i<exitCount; i++) {
            auto& it = exitInfoIterators[i];
            while((it != exitOrientedReadInfos[i].end()) and (it->orientedReadId < orientedReadId)) {
                ++it;
            }
            if((it == exitOrientedReadInfos[i].end()) or (it->orientedReadId != orientedReadId)) {
                commonOrientedReadInfo.exitStepCount[i] = 0;
            } else {
                commonOrientedReadInfo.exitStepCount[i] = it->stepCount;
                commonOrientedReadInfo.minPositionInJourneyOnExits =
                    min(commonOrientedReadInfo.minPositionInJourneyOnExits, it->positionInJourney);
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
    html << std::setprecision(2) << std::defaultfloat <<
        "<h4>Common oriented reads</h4>"
        "These are the oriented reads that appear in the representative region "
        "of at least one entrance and at least one exit. "
        "The ones that don't go backward contribute to the tangle matrix."
        "<p><table>"
        "<tr><th>Oriented<br>read<br>id";
    for(uint64_t i=0; i<entrances.size(); i++) {
        html << "<th>Entrance<br>" << assemblyGraph[entrances[i]].id;
    }
    for(uint64_t i=0; i<exits.size(); i++) {
        html << "<th>Exit<br>" << assemblyGraph[exits[i]].id;
    }
    html <<
        "<th>Max<br>position<br>in journey<br>(all entrances)"
        "<th>Min<br>position<br>in journey<br>(all exits)"
        "<th>Goes<br>backward?"
        "<th>Tangle<br>matrix";

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

        html <<
            "<td class=centered>" << commonOrientedReadInfo.maxPositionInJourneyOnEntrances <<
            "<td class=centered>" << commonOrientedReadInfo.minPositionInJourneyOnExits <<
            "<td class=centered>";
        if(commonOrientedReadInfo.goesBackward()) {
            html << "&#10003;";
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
        if(not commonOrientedReadInfo.goesBackward()) {
            for(uint64_t i=0; i<entranceCount; i++) {
                for(uint64_t j=0; j<exitCount; j++) {
                    tangleMatrix[i][j] += commonOrientedReadInfo.tangleMatrix[i][j];
                }
            }
        }
    }
}



void TangleMatrix1::writeTotalTangleMatrix(ostream& html) const
{
    html << "<h4>Total tangle matrix</h4><table><tr><th>";
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
