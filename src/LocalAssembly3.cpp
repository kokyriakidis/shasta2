#include "LocalAssembly3.hpp"
#include "Anchor.hpp"
#include "deduplicate.hpp"
#include "Markers.hpp"
#include "Reads.hpp"
using namespace shasta;



// This assembles between anchorIdA and anchorIdB
// of the given AnchorPair. It uses all the OrientedReadIds
// stored in the AnchorPair, and which all appear in both
// anchorIdA and anchorIdB. In addition, it uses OrientedReadIds
// stored in additionalOrientedReadIds that:
// - Are not also in the AnchorPair.
// - Appear in at least one of anchorIdA and anchorIdB.
LocalAssembly3::LocalAssembly3(
    const Anchors& anchors,
    ostream& html,
    bool debug,
    const AnchorPair& anchorPair,
    const vector<OrientedReadId>& additionalOrientedReadIds) :
    anchorIdA(anchorPair.anchorIdA),
    anchorIdB(anchorPair.anchorIdB)
{
    SHASTA_ASSERT(std::ranges::is_sorted(additionalOrientedReadIds));
    if(html) {
        writeInput(html, debug, anchorPair, additionalOrientedReadIds);
    }

    gatherOrientedReads(anchors,anchorPair, additionalOrientedReadIds);
    if(html) {
        writeOrientedReads(anchors, html);
    }
}



void LocalAssembly3::writeInput(
    ostream& html,
    bool debug,
    const AnchorPair& anchorPair,
    const vector<OrientedReadId>& additionalOrientedReadIds) const
{
    html <<
        "<table>"
        "<tr><th class=left>Left anchor<td class=centered>" << anchorIdToString(anchorIdA) <<
        "<tr><th class=left>Right anchor<td class=centered>" << anchorIdToString(anchorIdB) <<
        "</table>";

    if(debug) {
        html << "<h3>OrientedReadIds in the AnchorPair</h3><table>";
        for(const OrientedReadId orientedReadId: anchorPair.orientedReadIds) {
            html << "<tr><td class=centered>" << orientedReadId;
        }
        html << "</table>";

        html << "<h3>Additional OrientedReadIds</h3><table>";
        for(const OrientedReadId orientedReadId: additionalOrientedReadIds) {
            html << "<tr><td class=centered>" << orientedReadId;
        }
        html << "</table>";
    }

}



// Gather all the OrientedReadIds that we can consider for use in this assembly.
void LocalAssembly3::gatherOrientedReads(
    const Anchors& anchors,
    const AnchorPair& anchorPair,
    const vector<OrientedReadId>& additionalOrientedReadIds)
    {

    // Merge together the OrientedReadIds in the AnchorPair and
    // the additionalOrientedReadIds.
    vector<OrientedReadId> inputOrientedReadIds = anchorPair.orientedReadIds;
    std::ranges::copy(additionalOrientedReadIds, back_inserter(inputOrientedReadIds));
    deduplicate(inputOrientedReadIds);

    // Access the two Anchors.
    const Anchor anchorA = anchors[anchorIdA];
    const Anchor anchorB = anchors[anchorIdB];

    // To gather the OrientedReadIds that we can use for assembly,
    // use a joint loop over the inputOrientedReadIds and the oriented reads in
    // the two anchors.
    auto itA = anchorA.begin();
    const auto endA = anchorA.end();
    auto itB = anchorB.begin();
    const auto endB = anchorB.end();
    for(auto it=inputOrientedReadIds.begin(); it!=inputOrientedReadIds.end(); ++it) {
        const OrientedReadId orientedReadId = *it;

        while((itA !=endA) and (itA->orientedReadId < orientedReadId)) {
            ++itA;
        }
        while((itB != endB) and (itB->orientedReadId < orientedReadId)) {
            ++itB;
        }
        const bool isOnLeftAnchor = (itA != endA) and (itA->orientedReadId == orientedReadId);
        const bool isOnRightAnchor = (itB != endB) and (itB->orientedReadId == orientedReadId);

        if(isOnLeftAnchor or isOnRightAnchor) {
            const auto orientedReadMarkers = anchors.markers[orientedReadId.getValue()];

            OrientedReadInfo orientedReadInfo;
            orientedReadInfo.orientedReadId = orientedReadId;
            orientedReadInfo.isOnLeftAnchor = isOnLeftAnchor;
            orientedReadInfo.isOnRightAnchor = isOnRightAnchor;
            if(isOnLeftAnchor) {
                orientedReadInfo.leftOrdinal = itA->ordinal;
                orientedReadInfo.leftPosition = orientedReadMarkers[itA->ordinal].position;
            }
            if(isOnRightAnchor) {
                orientedReadInfo.rightOrdinal = itB->ordinal;
                orientedReadInfo.rightPosition = orientedReadMarkers[itB->ordinal].position;
            }
            orientedReadInfos.push_back(orientedReadInfo);
        }
    }
}



void LocalAssembly3::writeOrientedReads(
    const Anchors& anchors,
    ostream& html) const
{

    html << "<h3>OrientedReadIds used in this local assembly</h3>"
        "<table><tr>"
        "<th>Oriented<br>read id"
        "<th>Oriented<br>read<br>length"
        "<th>On left<br>anchor"
        "<th>On right<br>anchor"
        "<th>On both<br>anchors"
        "<th>Left<br>ordinal"
        "<th>Right<br>ordinal"
        "<th>Ordinal<br>offset"
        "<th>Left<br>position"
        "<th>Right<br>position"
        "<th>Position<br>offset";

    for(const OrientedReadInfo& orientedReadInfo: orientedReadInfos) {
        html <<
            "<tr><td class=centered>" << orientedReadInfo.orientedReadId <<
            "<td class=centered>" <<
            anchors.reads.getReadSequenceLength(orientedReadInfo.orientedReadId.getReadId());

        html << "<td class=centered>";
        if(orientedReadInfo.isOnLeftAnchor) {
            html << "&#10004;";
        }
        html << "<td class=centered>";
        if(orientedReadInfo.isOnRightAnchor) {
            html << "&#10004;";
        }
        html << "<td class=centered>";
        if(orientedReadInfo.isOnBothAnchors()) {
            html << "&#10004;";
        }
        html << "<td class=centered>";
        if(orientedReadInfo.isOnLeftAnchor) {
            html << orientedReadInfo.leftOrdinal;
        }
        html << "<td class=centered>";
        if(orientedReadInfo.isOnRightAnchor) {
            html << orientedReadInfo.rightOrdinal;
        }
        html << "<td class=centered>";
        if(orientedReadInfo.isOnBothAnchors()) {
            html << orientedReadInfo.ordinalOffset();
        }
        html << "<td class=centered>";
        if(orientedReadInfo.isOnLeftAnchor) {
            html << orientedReadInfo.leftPosition;
        }
        html << "<td class=centered>";
        if(orientedReadInfo.isOnRightAnchor) {
            html << orientedReadInfo.rightPosition;
        }
        html << "<td class=centered>";
        if(orientedReadInfo.isOnBothAnchors()) {
            html << orientedReadInfo.positionOffset();
        }
    }

}

