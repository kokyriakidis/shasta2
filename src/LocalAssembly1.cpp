#include "LocalAssembly1.hpp"
#include "abpoaWrapper.hpp"
#include "deduplicate.hpp"
#include "Markers.hpp"
#include "orderPairs.hpp"
#include "orderVectors.hpp"
#include "poastaWrapper.hpp"
#include "Reads.hpp"
using namespace shasta;

#include <chrono.hpp>


LocalAssembly1::LocalAssembly1(
    const Anchors& anchors,
    AnchorId anchorIdA,
    AnchorId anchorIdB,
    bool computeAlignment,
    uint64_t maxAbpoaLength,
    ostream& html) :
    anchors(anchors),
    html(html)
{
    gatherOrientedReads(anchorIdA, anchorIdB);
    if(html) {
        writeOrientedReads();
    }

    usePoasta = maximumLength() > maxAbpoaLength;
    if(usePoasta) {
        runPoasta();
    } else {
        runAbpoa(computeAlignment);
    }

    if(html) {
        html << "<p>MSA was computed using " << (usePoasta ? "poasta" : "abpoa");
        writeConsensus();
        if(computeAlignment) {
            writeAlignment();
        }
    }

}



uint64_t LocalAssembly1::maximumLength() const
{
    uint64_t length = 0;
    for(const OrientedRead& orientedRead: orientedReads) {
        length = max(length, uint64_t(orientedRead.sequenceLength()));
    }
    return length;
}



void LocalAssembly1::gatherOrientedReads(
    AnchorId anchorIdA,
    AnchorId anchorIdB)
{
    orientedReads.clear();
    const uint32_t kHalf = uint32_t(anchors.k / 2);
    const auto anchorA = anchors[anchorIdA];
    const auto anchorB = anchors[anchorIdB];
    const auto beginA = anchorA.begin();
    const auto beginB = anchorB.begin();
    const auto endA = anchorA.end();
    const auto endB = anchorB.end();
    auto itA = beginA;
    auto itB = beginB;
    while((itA != endA) and (itB != endB)) {

        if(itA->orientedReadId < itB->orientedReadId) {
            ++itA;
            continue;
        }

        if(itB->orientedReadId < itA->orientedReadId) {
            ++itB;
            continue;
        }

        // We found a common OrientedReadId.
        OrientedReadId orientedReadId = itA->orientedReadId;
        SHASTA_ASSERT(orientedReadId == itB->orientedReadId);

        // Only use it if the ordinal offset positive.
        if(itB->ordinal > itA->ordinal) {
            const auto orientedReadMarkers = anchors.markers[orientedReadId.getValue()];
            OrientedRead orientedRead;
            orientedRead.orientedReadId = orientedReadId;
            orientedRead.positionInJourneyA = itA->positionInJourney;
            orientedRead.positionInJourneyB = itB->positionInJourney;
            orientedRead.ordinalA = itA->ordinal;
            orientedRead.ordinalB = itB->ordinal;
            const Marker& markerA = orientedReadMarkers[orientedRead.ordinalA];
            const Marker& markerB = orientedReadMarkers[orientedRead.ordinalB];
            orientedRead.basePositionA = markerA.position + kHalf;
            orientedRead.basePositionB = markerB.position + kHalf;

            orientedReads.push_back(orientedRead);
        }

        ++itA;
        ++itB;

    }
    sort(orientedReads.begin(), orientedReads.end());


    // Now fill in the sequences.
    for(OrientedRead& orientedRead: orientedReads) {
        const OrientedReadId orientedReadId = orientedRead.orientedReadId;

        // Construct its sequence.
        for(uint32_t position=orientedRead.basePositionA; position!=orientedRead.basePositionB; position++) {
            orientedRead.sequence.push_back(anchors.reads.getOrientedReadBase(orientedReadId, position));
        }
    }
}



void LocalAssembly1::writeOrientedReads() const
{
    html <<
        "<h3>Oriented read sequences</h3>"
        "<p>This local assembly will use the following "<< orientedReads.size() << " oriented reads, "
        "ordered by decreasing sequence length."
        "<table>"
        "<tr>"
        "<th>OrientedReadId"
        "<th>PositionA<br>in journey<th>PositionB<br>in journey<th>Journey<br>offset"
        "<th>OrdinalA<th>OrdinalB<th>Ordinal<br>offset"
        "<th>PositionA<th>PositionB<th>Sequence<br>length"
        "<th>Sequence";

    for(const OrientedRead& orientedRead: orientedReads) {
        html <<
            "<tr>"
            "<td class=centered>" << orientedRead.orientedReadId <<
            "<td class=centered>" << orientedRead.positionInJourneyA <<
            "<td class=centered>" << orientedRead.positionInJourneyB <<
            "<td class=centered>" << orientedRead.positionInJourneyB - orientedRead.positionInJourneyA <<
            "<td class=centered>" << orientedRead.ordinalA <<
            "<td class=centered>" << orientedRead.ordinalB <<
            "<td class=centered>" << orientedRead.ordinalB - orientedRead.ordinalA <<
            "<td class=centered>" << orientedRead.basePositionA <<
            "<td class=centered>" << orientedRead.basePositionB <<
            "<td class=centered>" << orientedRead.basePositionB - orientedRead.basePositionA <<
            "<td style='font-family:monospace'>";
        copy(orientedRead.sequence.begin(), orientedRead.sequence.end(), ostream_iterator<Base>(html));
    }

    html << "</table>";
}


#if 0
void LocalAssembly1::gatherSequences()
{
    // Create a vector of indexes in the orientedReads vector
    // together with the corresponding sequence length.
    // Sort it by decreasing sequence length.
    vector< pair<uint64_t, uint64_t> > orientedReadTable;
    for(uint64_t orientedReadIndex=0; orientedReadIndex<orientedReads.size(); orientedReadIndex++) {
        const uint64_t sequenceLength = orientedReads[orientedReadIndex].sequenceLength();
        orientedReadTable.push_back(make_pair(orientedReadIndex, sequenceLength));
    }
    sort(orientedReadTable.begin(), orientedReadTable.end(), OrderPairsBySecondOnlyGreater<uint64_t, uint64_t>());



    // Process the OrientedReads in order of decreasing sequence length.
    vector<Base> orientedReadSequence;
    for(const auto& p: orientedReadTable) {
        const uint64_t orientedReadIndex = p.first;
        OrientedRead& orientedRead = orientedReads[orientedReadIndex];
        const OrientedReadId orientedReadId = orientedRead.orientedReadId;

        // Construct its sequence.
        orientedReadSequence.clear();
        for(uint32_t position=orientedRead.basePositionA; position!=orientedRead.basePositionB; position++) {
            orientedReadSequence.push_back(anchors.reads.getOrientedReadBase(orientedReadId, position));
        }
        SHASTA_ASSERT(orientedReadSequence.size() == p.second);

        // Look it up in the sequences we already stored.
        bool found = false;
        for(uint64_t sequenceId=0; sequenceId<sequences.size(); sequenceId++) {
            Sequence& sequence = sequences[sequenceId];
            if(sequence.sequence == orientedReadSequence) {
                sequence.orientedReadIndexes.push_back(orientedReadIndex);
                orientedRead.sequenceId = sequenceId;
                found = true;
                break;
            }
        }
        if(not found) {
            orientedRead.sequenceId = sequences.size();
            sequences.emplace_back(orientedReadSequence, orientedReadIndex);
        }
    }
}



void LocalAssembly1::writeSequences() const
{
    html << "<h3>Distinct sequences</h3>"
        "<table><tr><th>Sequence<br>id<th>Coverage<th>Length<th>Sequence";

    for(uint64_t sequenceId=0; sequenceId<sequences.size(); sequenceId++) {
        const Sequence& sequence = sequences[sequenceId];

        html << "<tr>"
            "<td class=centered>" << sequenceId <<
            "<td class=centered>" << sequence.orientedReadIndexes.size() <<
            "<td class=centered>" << sequence.sequence.size() <<
            "<td style='font-family:monospace'>";
        copy(sequence.sequence.begin(), sequence.sequence.end(), ostream_iterator<Base>(html));
    }

    html << "</table>";

}
#endif


void LocalAssembly1::runAbpoa(bool computeAlignment)
{
    vector< vector<Base> > inputSequences;

    for(const OrientedRead& orientedRead: orientedReads) {
        inputSequences.push_back(orientedRead.sequence);
    }

    const auto t0 = steady_clock::now();
    abpoa(inputSequences, consensus, alignment, alignedConsensus, computeAlignment);
    const auto t1 = steady_clock::now();
    cout << "abpoa ran in " << seconds(t1-t0) << " s." << endl;
}



void LocalAssembly1::runPoasta()
{
    vector< vector<Base> > inputSequences;

    for(const OrientedRead& orientedRead: orientedReads) {
        inputSequences.push_back(orientedRead.sequence);
    }

    const auto t0 = steady_clock::now();
    poasta(inputSequences, consensus, alignment, alignedConsensus);
    const auto t1 = steady_clock::now();
    cout << "poasta ran in " << seconds(t1-t0) << " s." << endl;
}



void LocalAssembly1::writeConsensus() const
{
    html <<
        "<h3>Consensus</h3>"
        "<table>"
        "<tr><th class=left>Consensus sequence length<td class=left>" << consensus.size() <<
        "<tr><th class=left>Consensus sequence"
        "<td style='font-family:monospace'>";

    for(uint64_t position=0; position<consensus.size(); position++) {
        const Base b = consensus[position].first;
        html << "<span title='" << position << "'>" << b << "</span>";
    }

    html <<
        "<tr><th class=left >Coverage"
        "<td style='font-family:monospace'>";

    std::map<char, uint64_t> coverageLegend;

    for(const auto& p: consensus) {
        const uint64_t coverage = p.second;
        const char c = (coverage < 10) ? char(coverage + '0') : char(coverage - 10 + 'A');
        coverageLegend.insert(make_pair(c, coverage));

        if(coverage < orientedReads.size()) {
            html << "<span style='background-color:Pink'>";
        }

        html << c;

        if(coverage < orientedReads.size()) {
            html << "</span>";
        }
    }

    html << "</table>";

    // Write the coverage legend.
    html << "<p><table><tr><th>Symbol<th>Coverage";
    for(const auto& p: coverageLegend) {
        html << "<tr><td class=centered>" << p.first << "<td class=centered>" << p.second;
    }
    html << "</table>";
}


void LocalAssembly1::writeAlignment() const
{
    html <<
        "<h3>Alignment</h3>"
        "<table>"
        "<tr><th class=left>OrientedReadId"
        "<th class=left>Sequence<br>length"
        "<th class=left>Aligned sequence";

    for(uint64_t i=0; i<alignment.size(); i++) {
        const vector<AlignedBase>& alignmentRow = alignment[i];

        html << "<tr><th>" << orientedReads[i].orientedReadId <<
            "<td class=centered>" << orientedReads[i].sequenceLength() <<
            "<td style='font-family:monospace;white-space: nowrap'>";

        for(uint64_t j=0; j<alignmentRow.size(); j++) {
            const AlignedBase b = alignmentRow[j];
            const bool isMatch = (b == alignedConsensus[j]);

            if(not isMatch) {
                html << "<span style='background-color:Pink'>";
            }
            html << b;
            if(not isMatch) {
                html << "</span>";
            }
        }

    }

    html << "<tr><th>Consensus<td>"
        "<td style='font-family:monospace;background-color:LightCyan;white-space:nowrap'>";

    uint64_t position = 0;
    for(uint64_t i=0; i<alignedConsensus.size(); i++) {
        const AlignedBase b = alignedConsensus[i];

        if(not b.isGap()) {
            html << "<span title='" << position << "'>";
        }

        html << b;

        if(not b.isGap()) {
            html << "</span>";
            ++position;
        }
    }

    html << "<tr><th>Consensus coverage<td>"
        "<td style='font-family:monospace;white-space:nowrap'>";

    position = 0;
    for(uint64_t i=0; i<alignedConsensus.size(); i++) {
        const AlignedBase b = alignedConsensus[i];

        if(b.isGap()) {
            html << "-";
        } else {
            const uint64_t coverage = consensus[position].second;
            const char c = (coverage < 10) ? char(coverage + '0') : char(coverage - 10 + 'A');

            if(coverage < orientedReads.size()) {
                html << "<span style='background-color:Pink'>";
            }

            html << c;

            if(coverage < orientedReads.size()) {
                html << "</span>";
            }

            ++position;
        }
    }

    html << "</table>";

}

