// Shasta.
#include "Anchor.hpp"
#include "color.hpp"
#include "deduplicate.hpp"
#include "ExternalAnchors.hpp"
#include "graphvizToHtml.hpp"
#include "hcsClustering.hpp"
#include "html.hpp"
#include "invalid.hpp"
#include "Journeys.hpp"
#include "MarkerInfo.hpp"
#include "MarkerKmers.hpp"
#include "Markers.hpp"
#include "orderPairs.hpp"
#include "orderVectors.hpp"
#include "performanceLog.hpp"
#include "Reads.hpp"
#include "runCommandWithTimeout.hpp"
#include "timestamp.hpp"
#include "tmpDirectory.hpp"
using namespace shasta2;

// Boost libraries.
#include <boost/dynamic_bitset.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>

// Standard library.
#include <cmath>
#include <queue>
#include <sstream>

// Explicit instantiation.
#include "MultithreadedObject.tpp"
template class MultithreadedObject<Anchors>;



// Constructor to accesses existing Anchors.
Anchors::Anchors(
    const string& baseName,
    const MappedMemoryOwner& mappedMemoryOwner,
    const Reads& reads,
    uint64_t k,
    const Markers& markers,
    const MarkerKmers& markerKmers,
    bool writeAccess) :
    MultithreadedObject<Anchors>(*this),
    MappedMemoryOwner(mappedMemoryOwner),
    baseName(baseName),
    reads(reads),
    k(k),
    kHalf(k/2),
    markers(markers),
    markerKmers(markerKmers)
{
    anchorMarkerInfos.accessExisting(largeDataName(baseName + "-AnchorMarkerInfos"), writeAccess);
    anchorInfos.accessExistingReadOnly(largeDataName(baseName + "-AnchorInfos"));
    kmerToAnchorTable.accessExistingReadOnly(largeDataName(baseName + "-KmerToAnchorTable"));
}



Anchor Anchors::operator[](AnchorId anchorId) const
{
    return anchorMarkerInfos[anchorId];
}



// This returns the sequence of the marker k-mer
// that this anchor was created from.
vector<Base> Anchors::anchorKmerSequence(AnchorId anchorId) const
{
    // Get the first AnchorMarkerInterval for this Anchor.
    const Anchor anchor = (*this)[anchorId];
    const AnchorMarkerInfo& firstMarkerInfo = anchor.front();

    // Get the OrientedReadId and the ordinals.
    const OrientedReadId orientedReadId = firstMarkerInfo.orientedReadId;
    const uint32_t ordinal = firstMarkerInfo.ordinal;

    // Access the markers of this OrientedReadId.
    const auto orientedReadMarkers = markers[orientedReadId.getValue()];

    const Marker& marker = orientedReadMarkers[ordinal];

    const uint32_t begin = marker.position;
    const uint32_t end = begin + uint32_t(k);

    vector<Base> sequence;
    for(uint32_t position=begin; position!=end; position++) {
        sequence.push_back(reads.getOrientedReadBase(orientedReadId, position));
    }

    return sequence;
}



Kmer Anchors::anchorKmer(AnchorId anchorId) const
{
    // Get the first AnchorMarkerInterval for this Anchor.
    const Anchor anchor = (*this)[anchorId];
    const AnchorMarkerInfo& firstMarkerInfo = anchor.front();

    // Get the OrientedReadId and the ordinal.
    const OrientedReadId orientedReadId = firstMarkerInfo.orientedReadId;
    const uint32_t ordinal = firstMarkerInfo.ordinal;

    return markers.getKmer(orientedReadId, ordinal);

}



uint64_t Anchors::size() const
{
    return anchorMarkerInfos.size();
}



void Anchors::check() const
{
    const Anchors& anchors = *this;

    for(AnchorId anchorId=0; anchorId<size(); anchorId++) {
        const Anchor& anchor = anchors[anchorId];
        anchor.check();
    }
}



void Anchor::check() const
{
    const Anchor& anchor = *this;

    for(uint64_t i=1; i<size(); i++) {
        SHASTA2_ASSERT(anchor[i-1].orientedReadId.getReadId() < anchor[i].orientedReadId.getReadId());
    }
}



// Return the number of common oriented reads between two Anchors,
// counting only oriented reads that have a greater ordinal on anchorId1
// than they have on anchorId0.
uint64_t Anchors::countCommon(
    AnchorId anchorId0,
    AnchorId anchorId1) const
{
    const Anchors& anchors = *this;
    const Anchor anchor0 = anchors[anchorId0];
    const Anchor anchor1 = anchors[anchorId1];

    auto it0 = anchor0.begin();
    auto it1 = anchor1.begin();

    const auto end0 = anchor0.end();
    const auto end1 = anchor1.end();

    uint64_t count = 0;
    while((it0 != end0) and (it1 != end1)) {
        const OrientedReadId orientedReadId0 = it0->orientedReadId;
        const OrientedReadId orientedReadId1 = it1->orientedReadId;
        if(orientedReadId0 < orientedReadId1) {
            ++it0;
        } else if(orientedReadId1 < orientedReadId0) {
            ++it1;
        } else {
            if(it0->ordinal < it1->ordinal) {
                ++count;
            }
            ++it0;
            ++it1;
        }
    }

    return count;
}



// Same as above, but also compute the average offset in bases.
uint64_t Anchors::countCommon(
    AnchorId anchorId0,
    AnchorId anchorId1,
    uint64_t& baseOffset) const
{
    const Anchors& anchors = *this;
    const Anchor anchor0 = anchors[anchorId0];
    const Anchor anchor1 = anchors[anchorId1];

    auto it0 = anchor0.begin();
    auto it1 = anchor1.begin();

    const auto end0 = anchor0.end();
    const auto end1 = anchor1.end();

    uint64_t count = 0;
    uint64_t sumBaseOffsets = 0;
    while((it0 != end0) and (it1 != end1)) {
        const OrientedReadId orientedReadId0 = it0->orientedReadId;
        const OrientedReadId orientedReadId1 = it1->orientedReadId;
        if(orientedReadId0 < orientedReadId1) {
            ++it0;
        } else if(orientedReadId1 < orientedReadId0) {
            ++it1;
        } else {

            // We found a common oriented read.
            const uint32_t ordinal0 = it0->ordinal;
            const uint32_t ordinal1 = it1->ordinal;
            if(ordinal0 < ordinal1) {
                ++count;

                const OrientedReadId orientedReadId = it0->orientedReadId;
                const auto orientedReadMarkers = markers[orientedReadId.getValue()];

                const uint64_t position0 = uint64_t(orientedReadMarkers[ordinal0].position);
                const uint64_t position1 = uint64_t(orientedReadMarkers[ordinal1].position);
                SHASTA2_ASSERT(position1 > position0);
                sumBaseOffsets += position1 - position0;
            }

            ++it0;
            ++it1;
        }
    }

    baseOffset = uint64_t(std::round(double(sumBaseOffsets) / double(count)));
    return count;
}



void Anchors::analyzeAnchorPair(
    AnchorId anchorIdA,
    AnchorId anchorIdB,
    AnchorPairInfo& info
    ) const
{
    const Anchors& anchors = *this;

    // Prepare for the joint loop over OrientedReadIds of the two Anchors.
    const Anchor anchorA = anchors[anchorIdA];
    const Anchor anchorB = anchors[anchorIdB];
    const auto beginA = anchorA.begin();
    const auto beginB = anchorB.begin();
    const auto endA = anchorA.end();
    const auto endB = anchorB.end();

    // Store the total number of OrientedReadIds on the two edges.
    info.totalA = endA - beginA;
    info.totalB = endB - beginB;


    // Joint loop over the MarkerIntervals of the two Anchors,
    // to count the common oreiented reads and compute average offsets.
    info.common = 0;
    int64_t sumMarkerOffsets = 0;
    int64_t sumBaseOffsets = 0;
    auto itA = beginA;
    auto itB = beginB;
    while(itA != endA and itB != endB) {

        if(itA->orientedReadId < itB->orientedReadId) {
            ++itA;
            continue;
        }

        if(itB->orientedReadId < itA->orientedReadId) {
            ++itB;
            continue;
        }

        // We found a common OrientedReadId.
        ++info.common;
        const OrientedReadId orientedReadId = itA->orientedReadId;
        const auto orientedReadMarkers = markers[orientedReadId.getValue()];

        // Compute the offset in markers.
        const uint32_t ordinalA = itA->ordinal;
        const uint32_t ordinalB = itB->ordinal;
        sumMarkerOffsets += int64_t(ordinalB) - int64_t(ordinalA);

        // Compute the offset in bases.
        const int64_t positionA = int64_t(orientedReadMarkers[ordinalA].position);
        const int64_t positionB = int64_t(orientedReadMarkers[ordinalB].position);
        sumBaseOffsets += positionB - positionA;

        // Continue the joint loop.
        ++itA;
        ++itB;

    }
    info.onlyA = info.totalA - info.common;
    info.onlyB = info.totalB - info.common;

    // If there are no common reads, this is all we can do.
    if(info.common == 0) {
        info.offsetInMarkers = invalid<int64_t>;
        info.offsetInBases = invalid<int64_t>;
        info.onlyAShort = invalid<uint64_t>;
        info.onlyBShort = invalid<uint64_t>;
        return;
    }

    // Compute the estimated offsets.
    info.offsetInMarkers = int64_t(std::round(double(sumMarkerOffsets) / double(info.common)));
    info.offsetInBases = int64_t(std::round(double(sumBaseOffsets) / double(info.common)));



    // Now do the joint loop again, and count the onlyA and onlyB oriented reads
    // that are too short to appear in the other edge.
    itA = beginA;
    itB = beginB;
    uint64_t onlyACheck = 0;
    uint64_t onlyBCheck = 0;
    info.onlyAShort = 0;
    info.onlyBShort = 0;
    while(true) {
        if(itA == endA and itB == endB) {
            break;
        }

        else if(itB == endB or ((itA!=endA) and (itA->orientedReadId < itB->orientedReadId))) {
            // This oriented read only appears in Anchor A.
            ++onlyACheck;
            const OrientedReadId orientedReadId = itA->orientedReadId;
            const auto orientedReadMarkers = markers[orientedReadId.getValue()];
            const int64_t lengthInBases = int64_t(reads.getReadSequenceLength(orientedReadId.getReadId()));

            const uint32_t ordinalA = itA->ordinal;
            const int64_t positionA = int64_t(orientedReadMarkers[ordinalA].position);

            // Find the hypothetical positions of anchor B, assuming the estimated base offset.
            const int64_t positionB = positionA + info.offsetInBases;

            // If this ends up outside the read, this counts as onlyAShort.
            if(positionB < 0 or positionB >= lengthInBases) {
                ++info.onlyAShort;
            }

            ++itA;
            continue;
        }

        else if(itA == endA or ((itB!=endB) and (itB->orientedReadId < itA->orientedReadId))) {
            // This oriented read only appears in Anchor B.
            ++onlyBCheck;
            const OrientedReadId orientedReadId = itB->orientedReadId;
            const auto orientedReadMarkers = markers[orientedReadId.getValue()];
            const int64_t lengthInBases = int64_t(reads.getReadSequenceLength(orientedReadId.getReadId()));

            // Get the positions of edge B in this oriented read.
            const uint32_t ordinalB = itB->ordinal;
            const int64_t positionB = int64_t(orientedReadMarkers[ordinalB].position);

            // Find the hypothetical positions of anchor A, assuming the estimated base offset.
            const int64_t positionA = positionB - info.offsetInBases;

            // If this ends up outside the read, this counts as onlyBShort.
            if(positionA < 0 or positionA >= lengthInBases) {
                ++info.onlyBShort;
            }

            ++itB;
            continue;
        }

        else {
            // This oriented read appears in both edges. In this loop, we
            // don't need to do anything.
            ++itA;
            ++itB;
        }
    }
    SHASTA2_ASSERT(onlyACheck == info.onlyA);
    SHASTA2_ASSERT(onlyBCheck == info.onlyB);
}



void Anchors::writeHtml(
    AnchorId anchorIdA,
    AnchorId anchorIdB,
    AnchorPairInfo& info,
    const Journeys& journeys,
    ostream& html) const
{
    const Anchors& anchors = *this;

    // Begin the summary table.
    html <<
        "<table>"
        "<tr><th><th>On<br>anchor A<th>On<br>anchor B";

    // Total.
    html <<
        "<tr><th class=left>Total ";
    writeInformationIcon(html, "The total number of oriented reads on each of the two anchors.");
    html << "<td class=centered>" << info.totalA << "<td class=centered>" << info.totalB;

    // Common.
    html << "<tr><th class=left>Common ";
    writeInformationIcon(html, "The number of common oriented reads between the two anchors.");
    html <<
        "<td class=centered colspan = 2>" << info.common;

    // Only.
    html <<
        "<tr><th class=left>Only ";
    writeInformationIcon(html, "The number of oriented reads that appear in one anchor but not the other.");
    html <<
        "<td class=centered>" << info.onlyA << "<td class=centered>" << info.onlyB;

    // The rest of the summary table can only be written if there are common reads.
    if(info.common > 0) {

        // Only, short.
        html <<
            "<tr><th class=left>Only, short ";
        writeInformationIcon(html, "The number of oriented reads that appear in one anchor only "
            " and are too short to appear on the other anchor, based on the estimated base offset.");
        html <<
            "<td class=centered>" << info.onlyAShort << "<td class=centered>" << info.onlyBShort;

        // Only, missing.
        html <<
            "<tr><th class=left>Only, missing ";
        writeInformationIcon(html, "The number of oriented reads that appear in one anchor only "
            " and are not too short to appear on the other anchor, based on the estimated base offset.");
        html <<
            "<td class=centered>" << info.onlyA - info.onlyAShort << "<td class=centered>" << info.onlyB - info.onlyBShort;
    }

    // End the summary table.
    html << "</table>";

    // Only write out the rest if there are common reads.
    if(info.common == 0) {
        return;
    }

    // Write the table with Jaccard similarities and estimated offsets.
    using std::fixed;
    using std::setprecision;
    html <<
        "<br><table>"
        "<tr><th class=left>Jaccard similarity<td class=centered>" <<
        fixed << setprecision(2) << info.jaccard() <<
        "<tr><th class=left>Corrected Jaccard similarity<td class=centered>" <<
        fixed << setprecision(2) << info.correctedJaccard() <<
        "<tr><th class=left>Estimated offset in markers<td class=centered>" << info.offsetInMarkers <<
        "<tr><th class=left>Estimated offset in bases<td class=centered>" << info.offsetInBases <<
        "</table>";



    // Write the details table.
    html <<
        "<br>In the following table, positions in red are hypothetical, based on the above "
        "estimated base offset."
        "<p><table>";

    // Header row.
    html <<
        "<tr>"
        "<th class=centered rowspan=2>Oriented<br>read id"
        "<th class=centered colspan=3>Length"
        "<th colspan=3>Anchor A"
        "<th colspan=3>Anchor B"
        "<th colspan=3>Offset"
        "<th rowspan=2>Classification"
        "<tr>"
        "<th>Bases"
        "<th>Markers"
        "<th>Anchors"
        "<th>Base<br>Position"
        "<th>Marker<br>ordinal"
        "<th>Position<br>in journey"
        "<th>Base<br>Position"
        "<th>Marker<br>ordinal"
        "<th>Position<br>in journey"
        "<th>Base<br>Position"
        "<th>Marker<br>ordinal"
        "<th>Position<br>in journey";

    // Prepare for the joint loop over OrientedReadIds of the two anchors.
    const auto markerIntervalsA = anchors[anchorIdA];
    const auto markerIntervalsB = anchors[anchorIdB];
    const auto beginA = markerIntervalsA.begin();
    const auto beginB = markerIntervalsB.begin();
    const auto endA = markerIntervalsA.end();
    const auto endB = markerIntervalsB.end();

    // Joint loop over the AnchorMarkerIntervals of the two Anchors.
    auto itA = beginA;
    auto itB = beginB;
    while(true) {
        if(itA == endA and itB == endB) {
            break;
        }

        else if(itB == endB or ((itA!=endA) and (itA->orientedReadId < itB->orientedReadId))) {
            // This oriented read only appears in Anchor A.
            const OrientedReadId orientedReadId = itA->orientedReadId;
            const auto orientedReadMarkers = markers[orientedReadId.getValue()];
            const int64_t lengthInBases = int64_t(reads.getReadSequenceLength(orientedReadId.getReadId()));
            const auto journey = journeys[orientedReadId];

            // Get the positions of Anchor A in this oriented read.
            const uint32_t ordinalA = itA->ordinal;
            const int64_t positionA = int64_t(orientedReadMarkers[ordinalA].position);

            // Find the hypothetical positions of Anchor B, assuming the estimated base offset.
            const int64_t positionB = positionA + info.offsetInBases;
            const bool isShort = positionB<0 or positionB >= lengthInBases;

            html <<
                "<tr><td class=centered>"
                "<a href='exploreRead?readId=" << orientedReadId.getReadId() <<
                "&strand=" << orientedReadId.getStrand() << "'>" << orientedReadId << "</a>"
                "<td class=centered>" << lengthInBases <<
                "<td class=centered>" << orientedReadMarkers.size() <<
                "<td class=centered>" << journey.size() <<
                "<td class=centered>" << positionA <<
                "<td class=centered>" << ordinalA <<
                "<td class=centered>" << itA->positionInJourney <<
                "<td class=centered style='color:Red'>" << positionB <<
                "<td>"
                "<td class=centered style='color:Red'>" << "<td><td><td>"
                "<td class=centered>OnlyA, " << (isShort ? "short" : "missing");

            ++itA;
            continue;
        }

        else if(itA == endA or ((itB!=endB) and (itB->orientedReadId < itA->orientedReadId))) {
            // This oriented read only appears in Anchor B.
            const OrientedReadId orientedReadId = itB->orientedReadId;
            const auto orientedReadMarkers = markers[orientedReadId.getValue()];
            const int64_t lengthInBases = int64_t(reads.getReadSequenceLength(orientedReadId.getReadId()));
            const auto journey = journeys[orientedReadId];

            // Get the positions of Anchor B in this oriented read.
            const uint32_t ordinalB = itB->ordinal;
            const int64_t positionB = int64_t(orientedReadMarkers[ordinalB].position);

            // Find the hypothetical positions of edge A, assuming the estimated base offset.
            const int64_t positionA = positionB - info.offsetInBases;
            const bool isShort = positionA<0 or positionA >= lengthInBases;

            html <<
                "<tr><td class=centered>"
                "<a href='exploreRead?readId=" << orientedReadId.getReadId() <<
                "&strand=" << orientedReadId.getStrand() << "'>" << orientedReadId << "</a>"
                "<td class=centered>" << lengthInBases <<
                "<td class=centered>" << orientedReadMarkers.size() <<
                "<td class=centered>" << journey.size() <<
                "<td class=centered style='color:Red'>" << positionA <<
                "<td><td>"
                "<td class=centered>" << positionB <<
                "<td class=centered>" << ordinalB <<
                "<td class=centered>" << itB->positionInJourney <<
                "<td class=centered>" << "<td><td>"
                "<td class=centered>OnlyB, " << (isShort ? "short" : "missing");

            ++itB;
            continue;
        }

        else {
            // This oriented read appears in both Anchors.
            const OrientedReadId orientedReadId = itA->orientedReadId;
            const auto orientedReadMarkers = markers[orientedReadId.getValue()];
            const int64_t lengthInBases = int64_t(reads.getReadSequenceLength(orientedReadId.getReadId()));
            const auto journey = journeys[orientedReadId];

            // Get the positions of Anchor A in this oriented read.
            const uint32_t ordinalA = itA->ordinal;
            const int64_t positionA = int64_t(orientedReadMarkers[ordinalA].position);

            // Get the positions of Anchor B in this oriented read.
            const uint32_t ordinalB = itB->ordinal;
            const int64_t positionB = int64_t(orientedReadMarkers[ordinalB].position);

            // Compute estimated offsets.
            const int64_t ordinalOffset = uint64_t(ordinalB) - uint64_t(ordinalA);
            const int64_t baseOffset = positionB - positionA;

            html <<
                "<tr><td class=centered>"
                "<a href='exploreRead?readId=" << orientedReadId.getReadId() <<
                "&strand=" << orientedReadId.getStrand() << "'>" << orientedReadId << "</a>"
                "<td class=centered>" << lengthInBases <<
                "<td class=centered>" << orientedReadMarkers.size() <<
                "<td class=centered>" << journey.size() <<
                "<td class=centered>" << positionA <<
                "<td class=centered>" << ordinalA <<
                "<td class=centered>" << itA->positionInJourney <<
                "<td class=centered>" << positionB <<
                "<td class=centered>" << ordinalB <<
                "<td class=centered>" << itB->positionInJourney <<
                "<td class=centered>" << baseOffset <<
                "<td class=centered>" << ordinalOffset <<
                "<td class=centered>" << int64_t(itB->positionInJourney) - int64_t(itA->positionInJourney) <<
                "<td class=centered>Common";

            ++itA;
            ++itB;
        }
    }

    // Finish the details table.
    html << "</table>";

}



// Anchors are numbered such that each pair of reverse complemented
// AnchorIds are numbered (n, n+1), where n is even, n = 2*m.
// We represent an AnchorId as a string as follows:
// - AnchorId n is represented as m+
// - AnchorId n+1 is represented as m-
// For example, the reverse complemented pair (150, 151) is represented as (75+, 75-).

string shasta2::anchorIdToString(AnchorId n)
{
    std::ostringstream s;

    const AnchorId m = (n >> 1);
    s << m;

    if(n & 1) {
        s << "-";
    } else {
        s << "+";
    }

    return s.str();
}



AnchorId shasta2::anchorIdFromString(const string& s)
{

    if(s.size() < 2) {
        return invalid<AnchorId>;
    }

    const char cLast = s.back();

    uint64_t lastBit;
    if(cLast == '+') {
        lastBit = 0;
    } else if(cLast == '-') {
        lastBit = 1;
    } else {
        return invalid<AnchorId>;
    }

    const uint64_t m = std::stoul(s.substr(0, s.size()-1));

    return 2* m + lastBit;
}



// For a given AnchorId, follow the read journeys forward by one step.
// Return a vector of the AnchorIds reached in this way.
// The count vector is the number of oriented reads each of the AnchorIds.
void Anchors::findChildren(
    const Journeys& journeys,
    AnchorId anchorId,
    vector<AnchorId>& children,
    vector<uint64_t>& count,
    uint64_t minCoverage) const
{
    children.clear();
    for(const auto& markerInfo: anchorMarkerInfos[anchorId]) {
        const OrientedReadId orientedReadId = markerInfo.orientedReadId;
        const auto journey = journeys[orientedReadId];
        const uint64_t position = markerInfo.positionInJourney;
        const uint64_t nextPosition = position + 1;
        if(nextPosition < journey.size()) {
            const AnchorId nextAnchorId = journey[nextPosition];
            children.push_back(nextAnchorId);
        }
    }

    deduplicateAndCountWithThreshold(children, count, minCoverage);
}



// For a given AnchorId, follow the read journeys backward by one step.
// Return a vector of the AnchorIds reached in this way.
// The count vector is the number of oriented reads each of the AnchorIds.
void Anchors::findParents(
    const Journeys& journeys,
    AnchorId anchorId,
    vector<AnchorId>& parents,
    vector<uint64_t>& count,
    uint64_t minCoverage) const
{
    parents.clear();
    for(const auto& markerInfo: anchorMarkerInfos[anchorId]) {
        const OrientedReadId orientedReadId = markerInfo.orientedReadId;
        const auto journey = journeys[orientedReadId];
        const uint64_t position = markerInfo.positionInJourney;
        if(position > 0) {
            const uint64_t previousPosition = position - 1;
            const AnchorId previousAnchorId = journey[previousPosition];
            parents.push_back(previousAnchorId);
        }
    }

    deduplicateAndCountWithThreshold(parents, count, minCoverage);
}



// Get the ordinal for the AnchorMarkerInfo corresponding to a
// given AnchorId and OrientedReadId.
// This asserts if the given AnchorId does not contain an AnchorMarkerInfo
// for the requested OrientedReadId.
uint32_t Anchors::getOrdinal(AnchorId anchorId, OrientedReadId orientedReadId) const
{
    for(const auto& markerInfo: anchorMarkerInfos[anchorId]) {
        if(markerInfo.orientedReadId == orientedReadId) {
            return markerInfo.ordinal;
        }
    }

    SHASTA2_ASSERT(0);
}



// Get the positionInJourney for the AnchorMarkerInfo corresponding to a
// given AnchorId and OrientedReadId.
// This asserts if the given AnchorId does not contain an AnchorMarkerInfo
// for the requested OrientedReadId.
uint32_t Anchors::getPositionInJourney(AnchorId anchorId, OrientedReadId orientedReadId) const
{
    for(const auto& markerInfo: anchorMarkerInfos[anchorId]) {
        if(markerInfo.orientedReadId == orientedReadId) {
            return markerInfo.positionInJourney;
        }
    }

    SHASTA2_ASSERT(0);
}


// Get the AnchorMarkerInfo corresponding to a given AnchorId and OrientedReadId.
// This asserts if the given AnchorId does not contain an AnchorMarkerInfo
// for the requested OrientedReadId.
const AnchorMarkerInfo& Anchors::getAnchorMarkerInfo(AnchorId anchorId, OrientedReadId orientedReadId) const
{
    for(const auto& markerInfo: anchorMarkerInfos[anchorId]) {
        if(markerInfo.orientedReadId == orientedReadId) {
            return markerInfo;
        }
    }

    SHASTA2_ASSERT(0);
}


// Find out if the given AnchorId contains the specified OrientedReadId.
bool Anchors::anchorContains(AnchorId anchorId, OrientedReadId orientedReadId) const
{
    for(const auto& markerInfo: anchorMarkerInfos[anchorId]) {
        if(markerInfo.orientedReadId == orientedReadId) {
            return true;
        }
    }

    return false;
}



void Anchors::writeCoverageHistogram() const
{
    vector<uint64_t> histogram;
    for(AnchorId anchorId=0; anchorId<anchorMarkerInfos.size(); anchorId++) {
        const uint64_t coverage = anchorMarkerInfos.size(anchorId);
        if(coverage >= histogram.size()) {
            histogram.resize(coverage + 1, 0);
        }
        ++histogram[coverage];
    }

    ofstream csv("AnchorCoverageHistogram.csv");
    csv << "Coverage,Frequency\n";
    for(uint64_t coverage=0; coverage<histogram.size(); coverage++) {
        csv << coverage << "," << histogram[coverage] << "\n";
    }
}



// Constructor to create Anchors from MarkerKmers.
Anchors::Anchors(
    const string& baseName,
    const MappedMemoryOwner& mappedMemoryOwner,
    const Reads& reads,
    uint64_t k,
    const Markers& markers,
    const MarkerKmers& markerKmers,
    uint64_t minAnchorCoverage,
    uint64_t maxAnchorCoverage,
    const vector<uint64_t>& maxAnchorRepeatLength,
    uint64_t threadCount) :
    MultithreadedObject<Anchors>(*this),
    MappedMemoryOwner(mappedMemoryOwner),
    baseName(baseName),
    reads(reads),
    k(k),
    kHalf(k/2),
    markers(markers),
    markerKmers(markerKmers)
{

    performanceLog << timestamp << "Anchor creation begins." << endl;

    // Adjust the numbers of threads, if necessary.
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    // Store arguments so all threads can see them.
    ConstructData& data = constructData;
    data.minAnchorCoverage = minAnchorCoverage;
    data.maxAnchorCoverage = maxAnchorCoverage;
    data.maxAnchorRepeatLength = maxAnchorRepeatLength;

    // During multithreaded pass 1 we loop over all marker k-mers
    // and for each one we find out if it can be used to generate
    // a pair of anchors or not. If it can be used,
    // we also fill in the coverage - that is,
    // the number of usable MarkerInfos that will go in each of the
    // two anchors.
    const uint64_t markerKmerCount = markerKmers.size();
    data.coverage.createNew(largeDataName("tmp-kmerToAnchorInfos"), largeDataPageSize);
    data.coverage.resize(markerKmerCount);
    const uint64_t batchSize = 1000;
    setupLoadBalancing(markerKmerCount, batchSize);
    runThreads(&Anchors::constructThreadFunctionPass1, threadCount);



    // Assign AnchorIds to marker k-mers and allocate space for
    // each anchor.
    anchorMarkerInfos.createNew(
            largeDataName(baseName + "-AnchorMarkerInfos"),
            largeDataPageSize);
    anchorInfos.createNew(largeDataName(baseName + "-AnchorInfos"), largeDataPageSize);
    kmerToAnchorTable.createNew(largeDataName(baseName + "-KmerToAnchorTable"), largeDataPageSize);
    kmerToAnchorTable.resize(markerKmerCount);
    AnchorId anchorId = 0;
    for(uint64_t kmerIndex=0; kmerIndex<markerKmerCount; kmerIndex++) {
        const uint64_t coverage = data.coverage[kmerIndex];
        if(coverage == 0) {
            // This k-mer does not generate any anchors.
            kmerToAnchorTable[kmerIndex] = invalid<AnchorId>;
        } else {
            // This k-mer generates two anchors.
            anchorMarkerInfos.appendVector(coverage);
            anchorMarkerInfos.appendVector(coverage);
            kmerToAnchorTable[kmerIndex] = anchorId;
            anchorInfos.push_back(AnchorInfo(kmerIndex));
            anchorInfos.push_back(AnchorInfo(kmerIndex));
            anchorId += 2;
        }
    }
    data.coverage.remove();
    const uint64_t anchorCount = anchorMarkerInfos.size();
    SHASTA2_ASSERT(anchorId == anchorCount);
    SHASTA2_ASSERT(anchorInfos.size() == anchorCount);



    // In pass 2 we fill in the AnchorMarkerInfos for each anchor.
    SHASTA2_ASSERT((batchSize % 2) == 0);
    setupLoadBalancing(anchorCount, batchSize);
    runThreads(&Anchors::constructThreadFunctionPass2, threadCount);



    cout << "Number of anchors per strand: " << anchorCount / 2 << endl;
    performanceLog << timestamp << "Anchor creation from marker kmers ends." << endl;

}




void Anchors::constructThreadFunctionPass1(uint64_t /* threadId */)
{

    ConstructData& data = constructData;
    const uint64_t minAnchorCoverage = data.minAnchorCoverage;
    const uint64_t maxAnchorCoverage = data.maxAnchorCoverage;
    const vector<uint64_t> maxAnchorRepeatLength = data.maxAnchorRepeatLength;

    // Loop over batches of marker Kmers assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over marker k-mers assigned to this batch.
        for(uint64_t kmerIndex=begin; kmerIndex!=end; kmerIndex++) {

            SHASTA2_ASSERT(data.coverage[kmerIndex] == 0);

            // Get the MarkerInfos for this marker Kmer.
            const span<const MarkerInfo> markerInfos = markerKmers[kmerIndex];

            // If coverage is too low, don't generate an anchor.
            if(markerInfos.size() < minAnchorCoverage) {
                continue;
            }

            // Check for high coverage using all of the marker infos.
            if(markerInfos.size() > maxAnchorCoverage) {
                continue;
            }

            // Count the usable MarkerInfos.
            // These are the ones for which the ReadId is different from the ReadId
            // of the previous and next MarkerInfo.
            uint64_t  usableMarkerInfosCount = 0;
            for(uint64_t i=0; i<markerInfos.size(); i++) {
                const MarkerInfo& markerInfo = markerInfos[i];
                bool isUsable = true;

                // Check if same ReadId of previous MarkerInfo.
                if(i != 0) {
                    isUsable =
                        isUsable and
                        (markerInfo.orientedReadId.getReadId() != markerInfos[i-1].orientedReadId.getReadId());
                }

                // Check if same ReadId of next MarkerInfo.
                if(i != markerInfos.size() - 1) {
                    isUsable =
                        isUsable and
                        (markerInfo.orientedReadId.getReadId() != markerInfos[i+1].orientedReadId.getReadId());
                }

                if(isUsable) {
                    ++usableMarkerInfosCount;
                }
            }

            if(markerInfos.size() - usableMarkerInfosCount > 0) {
                continue;
            }

            // If coverage is too low, don't generate an anchor.
            if(usableMarkerInfosCount < minAnchorCoverage) {
                continue;
            }

            // Check for repeats.
            bool skipDueToRepeats = false;
            const Kmer kmer = markerKmers.getKmer(markerInfos.front());
            for(uint64_t i=0; i<maxAnchorRepeatLength.size(); i++) {
                const uint64_t period = i + 1;
                const uint64_t maxAllowedCopyNumber = maxAnchorRepeatLength[i];
                if(kmer.countExactRepeatCopies(period, k) > maxAllowedCopyNumber) {
                    skipDueToRepeats = true;
                    break;
                }
            }
            if(skipDueToRepeats) {
                continue;
            }

            // If getting here, we will generate a pair of Anchors corresponding to this Kmer.
            data.coverage[kmerIndex] = usableMarkerInfosCount;
        }
    }
}



void Anchors::constructThreadFunctionPass2(uint64_t /* threadId */)
{

    ConstructData& data = constructData;
    const uint64_t minAnchorCoverage = data.minAnchorCoverage;
    const uint64_t maxAnchorCoverage = data.maxAnchorCoverage;


    // A vector used below and defined here to reduce memory allocation activity.
    // It will contain the MarkerInfos for a marker Kmer, excluding
    // the ones for which the same ReadId appears more than once in the same Kmer.
    // There are the ones that will be used to generate anchors.
    vector<MarkerInfo> usableMarkerInfos;

    // Loop over batches of AnchorIds assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over marker k-mers assigned to this batch.
        for(AnchorId anchorId=0; anchorId!=end; anchorId+=2) {
            const uint64_t kmerIndex = anchorInfos[anchorId].kmerIndex;

            // Get the MarkerInfos for this marker Kmer.
            const span<const MarkerInfo> markerInfos = markerKmers[kmerIndex];

            // We already checked for high coverage during pass 1.
            SHASTA2_ASSERT(markerInfos.size() <= maxAnchorCoverage);

            // Gather the usable MarkerInfos.
            // These are the ones for which the ReadId is different from the ReadId
            // of the previous and next MarkerInfo.
            usableMarkerInfos.clear();
            for(uint64_t i=0; i<markerInfos.size(); i++) {
                const MarkerInfo& markerInfo = markerInfos[i];
                bool isUsable = true;

                // Check if same ReadId of previous MarkerInfo.
                if(i != 0) {
                    isUsable =
                        isUsable and
                        (markerInfo.orientedReadId.getReadId() != markerInfos[i-1].orientedReadId.getReadId());
                }

                // Check if same ReadId of next MarkerInfo.
                if(i != markerInfos.size() - 1) {
                    isUsable =
                        isUsable and
                        (markerInfo.orientedReadId.getReadId() != markerInfos[i+1].orientedReadId.getReadId());
                }

                if(isUsable) {
                    usableMarkerInfos.push_back(markerInfo);
                }
            }

            if(markerInfos.size() - usableMarkerInfos.size() > 0) {
                continue;
            }

            // We already checked for low coverage durign pass1.
            SHASTA2_ASSERT(usableMarkerInfos.size() >= minAnchorCoverage);

            // Fill in the AnchorMarkerInfos for this anchor.
            const auto& anchorMarkerInfos0 = anchorMarkerInfos[anchorId];
            SHASTA2_ASSERT(anchorMarkerInfos0.size() == usableMarkerInfos.size());
            copy(usableMarkerInfos.begin(), usableMarkerInfos.end(), anchorMarkerInfos0.begin());

            // Reverse complement the usableMarkerInfos, then
            // generate the second anchor in the pair.
            for(MarkerInfo& markerInfo: usableMarkerInfos) {
                markerInfo = markerInfo.reverseComplement(markers);
            }
            const auto& anchorMarkerInfos1 = anchorMarkerInfos[anchorId + 1];
            SHASTA2_ASSERT(anchorMarkerInfos1.size() == usableMarkerInfos.size());
            copy(usableMarkerInfos.begin(), usableMarkerInfos.end(), anchorMarkerInfos1.begin());
        }
    }
}



// Cluster oriented reads in an anchor pair using their journey
// portions between AnchorIdA and AnchorIdB.
// Output to html if it is open.
// Returns clusters in order of decreasing length.
// Each cluster contains indices in AnchoirPair::orientedReadIds
// of the OrientedReadIds that belong to that cluster.
void Anchors::clusterAnchorPairOrientedReads(
    const AnchorPair& anchorPair,
    const Journeys& journeys,
    double clusteringMinJaccard,
    vector< vector<uint64_t> >& clusters,
    ostream& html) const
{
    const uint64_t orientedReadCount = anchorPair.size();

    // Get positions of AnchorIdA and AnchorIdB on the oriented reads og this AnchorPair.
    vector< pair<AnchorPair::Positions, AnchorPair::Positions> > positions;
    anchorPair.get(*this, positions);

    // Gather the AnchorIds visited by each oriented read between AnchorIdA and AnchorIdB
    // of the given AnchorPair.
    if(html) {
        html << "<h3>Journey portions within this anchor pair</h3><p><table>";
    }
    vector<AnchorId> anchorIds;
    for(uint64_t i=0; i<orientedReadCount; i++) {
        const OrientedReadId orientedReadId = anchorPair.orientedReadIds[i];
        const auto& positionsAB = positions[i];

        const Journey journey = journeys[orientedReadId];

        const auto& positionsA = positionsAB.first;
        const auto& positionsB = positionsAB.second;
        const auto positionInJourneyA = positionsA.positionInJourney;
        const auto positionInJourneyB = positionsB.positionInJourney;

        if(html) {
            html << "<tr><th class=centered>" << orientedReadId;
        }
        for(auto position=positionInJourneyA+1; position<positionInJourneyB; position++) {
            const AnchorId anchorId = journey[position];
            anchorIds.push_back(anchorId);
            html << "<td class=centered>" << anchorIdToString(journey[position]);
        }
    }
    deduplicate(anchorIds);
    const uint64_t anchorCount = anchorIds.size();

    if(html) {
        html << "</table>"
            "<p>Found " << anchorCount << " distinct anchors.";
    }

    // Create a bit vector for each OrienteRead, with a bit for each AnchorId.
    // The bit is set if the OrientedRead visits that AnchorId.
    using BitVector = boost::dynamic_bitset<uint64_t>;
    vector<BitVector> bitVectors(orientedReadCount);
    for(uint64_t i=0; i<orientedReadCount; i++) {
        const OrientedReadId orientedReadId = anchorPair.orientedReadIds[i];
        const auto& positionsAB = positions[i];

        const Journey journey = journeys[orientedReadId];

        const auto& positionsA = positionsAB.first;
        const auto& positionsB = positionsAB.second;
        const auto positionInJourneyA = positionsA.positionInJourney;
        const auto positionInJourneyB = positionsB.positionInJourney;

        BitVector& bitVector = bitVectors[i];
        bitVector.resize(anchorCount);

        for(auto position=positionInJourneyA+1; position<positionInJourneyB; position++) {
            const AnchorId anchorId = journey[position];
            const auto it = std::lower_bound(anchorIds.begin(), anchorIds.end(), anchorId);
            SHASTA2_ASSERT(it != anchorIds.end());
            SHASTA2_ASSERT(*it == anchorId);
            const uint64_t bitPosition = it - anchorIds.begin();
            bitVector.set(bitPosition);
        }
    }

    // Write out the bit vectors.
    if(html) {
        html << "<h3>OrientedReadId/AnchorId matrix</h3><p><table>";
        for(uint64_t i=0; i<orientedReadCount; i++) {
            const OrientedReadId orientedReadId = anchorPair.orientedReadIds[i];
            const BitVector& bitVector = bitVectors[i];
            html << "<tr>";
            for(uint64_t j=0; j<anchorCount; j++) {
                const AnchorId anchorId = anchorIds[j];
                const bool bitValue = bitVector[j];
                html << "<td class=centered";
                if(bitValue) {
                    html << " style='background-color:Green'";
                }
                html <<
                    " title='" << orientedReadId << " " << anchorIdToString(anchorId) << "'"
                    ">" << int(bitValue);
            }
        }
        html << "</table>";
    }


    // Compute a matrix of Jaccard similarities and a matrix of Hamming distances.
    vector< vector<double> > jaccard(orientedReadCount, vector<double>(orientedReadCount));
    vector< vector<uint64_t> > hamming(orientedReadCount, vector<uint64_t>(orientedReadCount));
    for(uint64_t i0=0; i0<orientedReadCount; i0++) {
        jaccard[i0][i0] = 1.;
        hamming[i0][i0] = 0;
        const BitVector& bitVector0 = bitVectors[i0];
        const uint64_t n0 = bitVector0.count();
        for(uint64_t i1=i0+1; i1<orientedReadCount; i1++) {
            const BitVector& bitVector1 = bitVectors[i1];
            const uint64_t n1 = bitVector1.count();
            const uint64_t intersectionCount = (bitVector0 & bitVector1).count();
            const uint64_t unionCount = n0 + n1 - intersectionCount;
            const double j = (unionCount == 0) ? 0. : double(intersectionCount) / double(unionCount);
            jaccard[i0][i1] = j;
            jaccard[i1][i0] = j;
            const uint64_t h = (bitVector0 ^ bitVector1).count();
            hamming[i0][i1] = h;
            hamming[i1][i0] = h;
        }
    }



    // Write out the Jaccard matrix and the Hamming distance matrix.
    if(html) {
        html << "<h3>Jaccard similarity matrix</h3><p><table>";
        for(uint64_t i0=0; i0<orientedReadCount; i0++) {
            const OrientedReadId orientedReadId0 = anchorPair.orientedReadIds[i0];
            html << "<tr>";
            for(uint64_t i1=0; i1<orientedReadCount; i1++) {
                const OrientedReadId orientedReadId1 = anchorPair.orientedReadIds[i1];
                const double j = jaccard[i0][i1];
                const double H = j / 3.;
                const string color = hslToRgbString(H, 0.5, 0.5);
                html << "<td class=centered style='background-color:" << color <<
                    "' title='" << orientedReadId0 << " " << orientedReadId1 <<
                    "'>" << std::fixed << std::setprecision(2) << j;
            }
        }
        html << "</table>";


        html << "<h3>Hamming distance matrix</h3><p><table>";
        for(uint64_t i0=0; i0<orientedReadCount; i0++) {
            const OrientedReadId orientedReadId0 = anchorPair.orientedReadIds[i0];
            html << "<tr>";
            for(uint64_t i1=0; i1<orientedReadCount; i1++) {
                const OrientedReadId orientedReadId1 = anchorPair.orientedReadIds[i1];
                html << "<td class=centered title='" <<
                    orientedReadId0 << " " << orientedReadId1 <<
                    "'>" << hamming[i0][i1];
            }
        }
        html << "</table>";
    }


    // Create a similarity graph with a vertex for each oriented read.
    // Jaccard values not smaller than minJaccard generate an edge.
    class SimilarityGraphVertex {
    public:
        uint64_t clusterId = invalid<uint64_t>;
    };
    using SimilarityGraph = boost::adjacency_list<
        boost::vecS,
        boost::vecS,
        boost::undirectedS,
        SimilarityGraphVertex>;
    SimilarityGraph similarityGraph(orientedReadCount);
    for(uint64_t i0=0; i0<orientedReadCount; i0++) {
        for(uint64_t i1=i0+1; i1<orientedReadCount; i1++) {
            if(jaccard[i0][i1] >= clusteringMinJaccard) {
                add_edge(i0, i1, similarityGraph);
            }
        }
    }

    // Clustering of the similarity graph.
    hcsClustering(similarityGraph, clusters);
    sort(clusters.begin(), clusters.end(), OrderVectorsByDecreasingSize<uint64_t>());
    for(uint64_t clusterId=0; clusterId<clusters.size(); clusterId++) {
        for(const SimilarityGraph::vertex_descriptor v: clusters[clusterId]) {
            similarityGraph[v].clusterId = clusterId;
        }
    }



    // Display the similarity graph.
    if(html) {

        // Write it out in Graphviz format.
        const string uuid = to_string(boost::uuids::random_generator()());
        const string dotFileName = tmpDirectory() + uuid + ".dot";
        ofstream dot(dotFileName);
        dot << "graph SimilarityGraph {\n";
        for(uint64_t i=0; i<orientedReadCount; i++) {
            const OrientedReadId orientedReadId = anchorPair.orientedReadIds[i];
            const string color = hslToRgbString(double(similarityGraph[i].clusterId) / double(clusters.size()), 0.75, 0.6);
            dot << i << " [label=\"" << orientedReadId << "\\n" << similarityGraph[i].clusterId << "\""
                " fillcolor=\"" << color << "\"];\n";
        }
        BGL_FORALL_EDGES(e, similarityGraph, SimilarityGraph) {
            const uint64_t i0 = source(e, similarityGraph);
            const uint64_t i1 = target(e, similarityGraph);
            dot << i0 << "--" << i1 << ";\n";
        }
        dot << "}\n";
        dot.close();

        // Send it to html in svg format.
        const double timeout = 30.;
        const string options = "-Nshape=rectangle -Nstyle=filled -Goverlap=false -Gsplines=true -Gbgcolor=gray95";
        graphvizToHtml(dotFileName, "sfdp", timeout, options, html);
    }


}



// Use the kmerToAnchorTable table to get the AnchorId corresponding to a given Kmer.
// When using ExternalAnchors, this always returns invalid<AnchorId>.
// This is only used in the http server.
// It is not used in the standard assembly process.
AnchorId Anchors::getAnchorIdFromKmer(const Kmer& kmer) const
{
    const Kmer kmerRc = kmer.reverseComplement(k);

    if(kmer <= kmerRc) {

        // kmer is canonical.

        const uint64_t globalIndex = markerKmers.getGlobalIndex(kmer);
        const uint64_t anchorId = kmerToAnchorTable[globalIndex];
        if(anchorId == invalid<uint64_t>) {
            return invalid<uint64_t>;
        } else {
            return anchorId;
        }

    } else {

        // kmerRc is canonical.

        const uint64_t globalIndex = markerKmers.getGlobalIndex(kmerRc);
        const uint64_t anchorId = kmerToAnchorTable[globalIndex];
        if(anchorId == invalid<uint64_t>) {
            return invalid<uint64_t>;
        } else {
            return anchorId + 1;
        }

    }
}



// Constructor to read Anchors from ExternalAnchors.
Anchors::Anchors(
    const string& baseName,
    const MappedMemoryOwner& mappedMemoryOwner,
    const Reads& reads,
    uint64_t k,
    const Markers& markers,
    const MarkerKmers& markerKmers,
    const string& externalAnchorsName) :
    MultithreadedObject<Anchors>(*this),
    MappedMemoryOwner(mappedMemoryOwner),
    baseName(baseName),
    reads(reads),
    k(k),
    kHalf(k/2),
    markers(markers),
    markerKmers(markerKmers)
{

    // Access the ExternalAnchors.
    SHASTA2_ASSERT(not externalAnchorsName.empty());
    if(externalAnchorsName[0] != '/') {
        throw runtime_error("--external-anchors-name must specify an absolute path.");
    }
    cout << "Reading external anchors " << externalAnchorsName << endl;
    const ExternalAnchors externalAnchors(externalAnchorsName, ExternalAnchors::AccessExisting());
    cout << "Found " << externalAnchors.data.size() <<
        " external anchors with average coverage " <<
        externalAnchors.data.totalSize() / externalAnchors.data.size() << endl;

    // Initialize the binary data owned by Anchors.
    anchorMarkerInfos.createNew(
        largeDataName(baseName + "-AnchorMarkerInfos"),
        largeDataPageSize);
    anchorInfos.createNew(largeDataName(baseName + "-AnchorInfos"), largeDataPageSize);

    // Create the kmerToAnchorTable and fill it in with invalid<AnchorId>.
    kmerToAnchorTable.createNew(largeDataName(baseName + "-KmerToAnchorTable"), largeDataPageSize);
    kmerToAnchorTable.resize(markerKmers.size());
    fill(kmerToAnchorTable.begin(), kmerToAnchorTable.end(), invalid<AnchorId>);


    // Loop over external anchors.
    // Each external anchor generates a pair of Anchors.
    vector<AnchorMarkerInfo> markerInfos;
    for(uint64_t i=0; i<externalAnchors.data.size(); i++) {
        const span<const ExternalAnchors::OrientedRead> externalAnchor = externalAnchors.data[i];

        // Create the MarkerInfos for this external anchor.
        // Also check that the Kmers for all oriented reads are identical.
        Kmer kmer;
        markerInfos.clear();
        for(const ExternalAnchors::OrientedRead& orientedRead: externalAnchor) {
            const OrientedReadId orientedReadId = orientedRead.orientedReadId;
            const uint32_t position = orientedRead.position;
            const span<const Marker> orientedReadMarkers = markers[orientedReadId.getValue()];

            // Locate the marker at this position.
            Marker targetMarker;
            targetMarker.position = position;
            const auto it = std::lower_bound(orientedReadMarkers.begin(), orientedReadMarkers.end(), targetMarker);
            if((it==orientedReadMarkers.end()) or (it->position != position)) {
                std::ostringstream message;
                message << "Oriented read " << orientedReadId <<
                    " does not have a marker at position " << position;
                cout << message.str() << endl;
                cout << "Offending external anchor:" << endl;
                externalAnchors.write(cout, i, k, reads, markers);
                throw runtime_error(message.str());
            }
            const uint32_t ordinal = uint32_t(it - orientedReadMarkers.begin());

            // Check the Kmer.
            const Kmer orientedReadKmer = markers.getKmer(orientedReadId, ordinal);
            if(markerInfos.empty()) {
                kmer = orientedReadKmer;
            } else {
                if(orientedReadKmer != kmer) {
                    std::ostringstream message;
                    message << "Inconsistent kmer at oriented read " << orientedReadId <<
                        " position " << position << " ordinal " << ordinal << endl;
                    cout << message.str() << endl;
                    cout << "Offending external anchor:" << endl;
                    externalAnchors.write(cout, i, k, reads, markers);
                    throw runtime_error(message.str());
                }
            }

            // Store this MarkerInfo.
            markerInfos.emplace_back(orientedReadId, ordinal);
        }

        // Sort the MarkerInfos by OrientedReadId.
        sort(markerInfos.begin(), markerInfos.end());

        // Check that the ReadIds are all distinct.
        for(uint64_t i1=1; i1<markerInfos.size(); i1++) {
            const uint64_t i0 = i1 - 1;
            if(markerInfos[i0].orientedReadId.getReadId() >= markerInfos[i1].orientedReadId.getReadId()) {
                std::ostringstream message;
                message << "Duplicate ReadId " << markerInfos[i0].orientedReadId.getReadId() << endl;
                cout << "Offending external anchor:" << endl;
                externalAnchors.write(cout, i, k, reads, markers);
                throw runtime_error(message.str());
            }
        }

        // Generate the Anchor corresponding to this external anchor,
        // without reverse complementing.
        // We are not filling AnchorInfo::kmerIndex.
        anchorMarkerInfos.appendVector(markerInfos);
        anchorInfos.push_back(AnchorInfo());

        // Reverse complement the MarkerInfos, then generate
        // the reverse complemented Anchor.
        // We are not filling AnchorInfo::kmerIndex.
        for(AnchorMarkerInfo& markerInfo: markerInfos) {
            const uint32_t orientedReadMarkerCount = uint32_t(markers.size(markerInfo.orientedReadId.getValue()));
            markerInfo.orientedReadId.flipStrand();
            markerInfo.ordinal = orientedReadMarkerCount - 1 - markerInfo.ordinal;
        }
        anchorMarkerInfos.appendVector(markerInfos);
        anchorInfos.push_back(AnchorInfo());
    }

    cout << "Generated " << anchorMarkerInfos.size() << " anchors from " <<
        externalAnchors.data.size() << " external anchors." << endl;
}



// Constructor to create an empty Anchors object.
Anchors::Anchors(
    const string& baseName,
    const MappedMemoryOwner& mappedMemoryOwner,
    const Reads& reads,
    uint64_t k,
    const Markers& markers,
    const MarkerKmers& markerKmers) :
    MultithreadedObject<Anchors>(*this),
    MappedMemoryOwner(mappedMemoryOwner),
    baseName(baseName),
    reads(reads),
    k(k),
    kHalf(k/2),
    markers(markers),
    markerKmers(markerKmers)
{
    anchorMarkerInfos.createNew(largeDataName(baseName + "-AnchorMarkerInfos"), largeDataPageSize);
    anchorInfos.createNew(largeDataName(baseName + "-AnchorInfos"), largeDataPageSize);
    kmerToAnchorTable.createNew(largeDataName(baseName + "-KmerToAnchorTable"), largeDataPageSize);
}



// Constructor that makes a copy of a source Anchors object, but removing
// a specified set of AnchorMarkerInfos.
// The keep vector specifies which AnchorMarkerInfos should be kept.
// It must be of size that.anchorMarkerInfos.totalSize() and
// is indexed by the global position of the AnchorMarkerInfo
// in that.anchorMarkerInfos, that is, &anchorMarkerInfo-that.anchorMarkerInfos.begin().
Anchors::Anchors(
    const Anchors& sourceAnchors,
    const string& baseName,
    const vector<bool>& keep) :
    MultithreadedObject<Anchors>(*this),
    MappedMemoryOwner(sourceAnchors),
    baseName(baseName),
    reads(sourceAnchors.reads),
    k(sourceAnchors.k),
    kHalf(k/2),
    markers(sourceAnchors.markers),
    markerKmers(sourceAnchors.markerKmers)
{
    // Initialize MemoryMapped objects.
    anchorMarkerInfos.createNew(largeDataName(baseName + "-AnchorMarkerInfos"), largeDataPageSize);
    anchorInfos.createNew(largeDataName(baseName + "-AnchorInfos"), largeDataPageSize);
    kmerToAnchorTable.createNew(largeDataName(baseName + "-KmerToAnchorTable"), largeDataPageSize);

    // Initialize the kmerToAnchorTable.
    kmerToAnchorTable.resize(sourceAnchors.markerKmers.size());
    std::ranges::fill(kmerToAnchorTable, invalid<AnchorId>);

    // Loop over anchors in the sourceAnchors.
    // Each of them generates an Anchor,
    // with some AnchorMarkerInfos possibly excluded.
    for(AnchorId sourceAnchorId=0; sourceAnchorId<sourceAnchors.size(); sourceAnchorId++) {
        const Anchor sourceAnchor = sourceAnchors[sourceAnchorId];

        // Get the anchorId for the new anchor we are generating.
        // Since we are nto removing any anchors, this must be the same
        // as the sourceAnchorId.
        const AnchorId anchorId = size();
        SHASTA2_ASSERT(anchorId == sourceAnchorId);

        // Copy the AnchorInfo.
        const AnchorInfo& sourceAnchorInfo = sourceAnchors.anchorInfos[sourceAnchorId];
        anchorInfos.push_back(sourceAnchorInfo);

        // Update the kmerToAnchorTable.
        kmerToAnchorTable[sourceAnchorInfo.kmerIndex] = anchorId;

        // Generate the AnchorMarkerInfos for this new anchor.
        // They are the same as for sourceAnchor, but some
        // AnchorMarkerInfos may be removed.
        anchorMarkerInfos.appendVector();
        uint64_t newCoverage = 0;
        for(const AnchorMarkerInfo& sourceAnchorMarkerInfo: sourceAnchor) {
            if(keep[&sourceAnchorMarkerInfo - sourceAnchors.anchorMarkerInfos.begin()]) {
                anchorMarkerInfos.append(sourceAnchorMarkerInfo);
                ++newCoverage;
            }
        }
        if(newCoverage != sourceAnchor.size()) {
            cout << anchorIdToString(anchorId) << " old coverage " << sourceAnchor.size() <<
                ", new coverage " << newCoverage << endl;
        }
    }
}



void Anchors::remove()
{
    anchorMarkerInfos.remove();
    anchorInfos.remove();
    kmerToAnchorTable.remove();
}
