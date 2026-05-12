// Shasta.
#include "LocalAssembly6.hpp"
#include "Anchor.hpp"
#include "Base.hpp"
#include "Markers.hpp"
#include "Reads.hpp"
#include "theseusWrapper.hpp"
using namespace shasta2;

// Standard library.
#include <cmath>
#include <map>



LocalAssembly6::LocalAssembly6(
    const Anchors& anchors,
    AnchorId anchorIdA,
    AnchorId anchorIdB,
    ostream& html,
    bool debug,
    const vector<OrientedReadId>& orientedReadIds) :
    anchors(anchors),
    anchorIdA(anchorIdA),
    anchorIdB(anchorIdB),
    html(html),
    debug(debug)
{
    if(html) {
        html << "<h3>Local assembly between " <<
            anchorIdToString(anchorIdA) << " and " <<
            anchorIdToString(anchorIdB) << "</h3>";
    }

    gatherOrientedReads(orientedReadIds);
    estimateOffset();
    gatherOrientedReadsSequences();
    writeOrientedReads();
    assemble();
}



void LocalAssembly6::gatherOrientedReads(
    const vector<OrientedReadId>& orientedReadIds)
{
    const Anchor anchorA = anchors[anchorIdA];
    const Anchor anchorB = anchors[anchorIdB];
    const uint32_t kHalf = uint32_t(anchors.k / 2);

    if(html) {
        html << "<h4>" << orientedReadIds.size() << " input oriented reads</h4><table>";
        for(const OrientedReadId orientedReadId: orientedReadIds) {
            html << "<tr><td class=centered>" << orientedReadId;
        }
        html << "</table>";
    }



    // Joint loop over the input OrientedReadIds and
    // the OrientedReadIds of the two Anchors.
    auto itA = anchorA.begin();
    const auto endA = anchorA.end();
    auto itB = anchorB.begin();
    const auto endB = anchorB.end();
    for(const OrientedReadId orientedReadId: orientedReadIds) {

        // Check if this OrientedReadId is on the two anchors.
        while((itA != endA) and (itA->orientedReadId < orientedReadId)) {
            ++itA;
        }
        while((itB != endB) and (itB->orientedReadId < orientedReadId)) {
            ++itB;
        }
        const bool isOnA = (itA != endA) and (itA->orientedReadId == orientedReadId);
        const bool isOnB = (itB != endB) and (itB->orientedReadId == orientedReadId);

        // If on neither anchor, this OrientedReadId cannot be used.
        if(not (isOnA or isOnB)) {
            continue;
        }

        // We can use this OrientedReadId for assembly.
        // Get its markers.
        const auto markers = anchors.markers[orientedReadId.getValue()];

        // Create an OrientedReadInfo.
        OrientedReadInfo& info = orientedReadInfos.emplace_back();
        info.orientedReadId = orientedReadId;
        if(isOnA) {
            info.ordinalA = itA->ordinal;
            info.positionA = markers[info.ordinalA].position + kHalf;
        }
        if(isOnB) {
            info.ordinalB = itB->ordinal;
            info.positionB = markers[info.ordinalB].position + kHalf;
        }

    }
}


void LocalAssembly6::estimateOffset()
{
    uint64_t sum = 0;
    uint64_t n = 0;
    for(const OrientedReadInfo& orientedReadInfo: orientedReadInfos) {
        if(orientedReadInfo.isOnBothAnchors()) {
            sum += orientedReadInfo.positionOffsetAB();
            ++n;
        }
    }

    SHASTA2_ASSERT(n > 0);
    offset = uint32_t(std::round(double(sum) / double(n)));
}



void LocalAssembly6::gatherOrientedReadsSequences()
{
    // For reads fixed on one side only, we use a sequence length
    // equal to offset + aDrift * offset + bDrift.
    const uint32_t length = offset + uint32_t(std::round(aDrift * double(offset) + bDrift));


    // Fill in the beginPosition and endPosition of each OrientedReadInfo.
    // Those are the position ranges of the sequences that will be used for assembly.
    for(OrientedReadInfo& info: orientedReadInfos) {

        if(info.isOnBothAnchors()) {
            info.positionBegin = info.positionA;
            info.positionEnd  = info.positionB;
        } else if(info.isOnAnchorA()) {
            SHASTA2_ASSERT(not info.isOnAnchorB());
            const ReadId readId = info.orientedReadId.getReadId();
            const uint32_t readLength = uint32_t(anchors.reads.getReadSequenceLength(readId));
            info.positionBegin = info.positionA;
            info.positionEnd  = min(readLength, info.positionBegin + length);
        } else if(info.isOnAnchorB()) {
            SHASTA2_ASSERT(not info.isOnAnchorA());
            info.positionEnd = info.positionB;
            if(info.positionEnd > length) {
                info.positionBegin = info.positionEnd - length;
            } else {
                info.positionBegin  = 0;
            }
        } else {
            SHASTA2_ASSERT(0);
        }
    }



    // Now we can create the sequence tables.
    const Reads& reads = anchors.reads;
    vector<Base> sequence;
    for(OrientedReadInfo& info: orientedReadInfos) {
        const OrientedReadId orientedReadId = info.orientedReadId;

        // Create the sequence of this read (portion to be used for this local assembly).
        sequence.clear();
        for(uint32_t position=info.positionBegin; position!=info.positionEnd; position++) {
            sequence.push_back(reads.getOrientedReadBase(orientedReadId, position));
        }

        // Access the sequence table appropriate for this oriented read.
        vector<SequenceInfo>& sequenceTable =
            (info.isOnBothAnchors() ? fixedSequencesTable :
            (info.isOnAnchorA() ? leftFixedSequencesTable : rightFixedSequencesTable));

        // See if we already have this sequence in this table.
        bool found = false;
        for(SequenceInfo& sequenceInfo: sequenceTable) {
            if(*(sequenceInfo.sequencePointer) == sequence) {
                sequenceInfo.orientedReadIds.push_back(orientedReadId);
                found = true;
                break;
            }
        }
        if(not found) {
            sequenceTable.emplace_back(orientedReadId, sequence);
        }

    }

    // Sort the sequence tables by decreasing coverage.
    sort(fixedSequencesTable.begin(), fixedSequencesTable.end());
    sort(leftFixedSequencesTable.begin(), leftFixedSequencesTable.end());
    sort(rightFixedSequencesTable.begin(), rightFixedSequencesTable.end());

    // Assign ids to the sequences.
    std::map<OrientedReadId, uint64_t> sequenceMap; // Key=OrientedReadId. Value = SequenceId.
    uint64_t sequenceId = 0;
    for(SequenceInfo& sequenceInfo: fixedSequencesTable) {
        sequenceInfo.id = sequenceId++;
        for(const OrientedReadId orientedReadId: sequenceInfo.orientedReadIds) {
            sequenceMap.insert(make_pair(orientedReadId, sequenceInfo.id));
        }
    }
    for(SequenceInfo& sequenceInfo: leftFixedSequencesTable) {
        sequenceInfo.id = sequenceId++;
        for(const OrientedReadId orientedReadId: sequenceInfo.orientedReadIds) {
            sequenceMap.insert(make_pair(orientedReadId, sequenceInfo.id));
        }
    }
    for(SequenceInfo& sequenceInfo: rightFixedSequencesTable) {
        sequenceInfo.id = sequenceId++;
        for(const OrientedReadId orientedReadId: sequenceInfo.orientedReadIds) {
            sequenceMap.insert(make_pair(orientedReadId, sequenceInfo.id));
        }
    }

    // Now we can fill in the sequenceId fields of the oOrientedReadInfos.
    for(OrientedReadInfo& info: orientedReadInfos) {
        info.sequenceId = sequenceMap.at(info.orientedReadId);
    }
}



void LocalAssembly6::writeOrientedReads() const
{
    if(not html) {
        return;
    }

    html << "<h4>" << orientedReadInfos.size() << " usable oriented reads</h4>";

    html << "<table><tr>"
        "<th>Oriented<br>read id"
        "<th>On A"
        "<th>On B"
        "<th>OrdinalA"
        "<th>OrdinalB"
        "<th>Marker<br>count"
        "<th>Ordinal<br>offset"
        "<th>PositionA"
        "<th>PositionB"
        "<th>Length<br>(bases)"
        "<th>PositionAB<br>offset"
        "<th>Assembly<br>position<br>begin"
        "<th>Assembly<br>position<br>end"
        "<th>Assembly<br>position<br>offset"
        "<th>Sequence<br>id"
        ;


    uint64_t commonCount = 0;
    for(const OrientedReadInfo& info: orientedReadInfos) {
        if(info.isOnBothAnchors()) {
            ++commonCount;
        }
        html <<
            "<tr>"
            "<th class=centered>" << info.orientedReadId <<
            "<td class=centered>" << (info.isOnAnchorA() ? "&check;" : "") <<
            "<td class=centered>" << (info.isOnAnchorB() ? "&check;" : "") <<
            "<td class=centered>" << (info.isOnAnchorA() ? to_string(info.ordinalA) : "") <<
            "<td class=centered>" << (info.isOnAnchorB() ? to_string(info.ordinalB) : "") <<
            "<td class=centered>" << anchors.markers[info.orientedReadId.getValue()].size() <<
            "<td class=centered>" << (info.isOnBothAnchors() ? to_string(info.ordinalOffsetAB()) : "") <<
            "<td class=centered>" << (info.isOnAnchorA() ? to_string(info.positionA) : "") <<
            "<td class=centered>" << (info.isOnAnchorB() ? to_string(info.positionB) : "") <<
            "<td class=centered>" << anchors.reads.getReadSequenceLength(info.orientedReadId.getReadId()) <<
            "<td class=centered>" << (info.isOnBothAnchors() ? to_string(info.positionOffsetAB()) : "") <<
            "<td class=centered>" << info.positionBegin <<
            "<td class=centered>" << info.positionEnd <<
            "<td class=centered>" << info.positionOffsetForAssembly() <<
            "<td class=centered>" << info.sequenceId
            ;
    }

    html << "</table>";

    html << "<br>Estimated offset using " << commonCount <<
        " oriented reads common to the left and right anchors is " << offset << " bases.";



    html << "<h4>Oriented read sequences used for assembly</h4>"
        "<table><tr>"
        "<th>Sequence<br>id"
        "<th>Fixed<br>on A"
        "<th>Fixed<br>on B"
        "<th>Length"
        "<th>Coverage"
        "<th class=left>Sequence";
    writeSequenceTable(fixedSequencesTable, true, true);
    writeSequenceTable(leftFixedSequencesTable, true, false);
    writeSequenceTable(rightFixedSequencesTable, false, true);

    html << "</table>";

}



void LocalAssembly6::writeSequenceTable(
    const vector<SequenceInfo>& sequenceTable,
    bool fixedOnA,
    bool fixedOnB) const
{
    for(const SequenceInfo& sequenceInfo: sequenceTable) {
        html <<
            "<tr>"
            "<td class=centered>" << sequenceInfo.id <<
            "<td class=centered>" << (fixedOnA ? "&check;" : "") <<
            "<td class=centered>" << (fixedOnB ? "&check;" : "") <<
            "<td class=centered>" << sequenceInfo.sequencePointer->size() <<
            "<td class=centered>" << sequenceInfo.orientedReadIds.size() <<
            "<td class=left style='font-family:monospace'>";
        std::ranges::copy(*(sequenceInfo.sequencePointer), ostream_iterator<Base>(html));
    }
}



void LocalAssembly6::assemble()
{
    vector< pair<vector<Base>, uint64_t> > fixedSequences;
    vector< pair<vector<Base>, uint64_t> > leftFixedSequences;
    vector< pair<vector<Base>, uint64_t> > rightFixedSequences;

    for(const SequenceInfo& sequenceInfo: fixedSequencesTable) {
        fixedSequences.emplace_back(*(sequenceInfo.sequencePointer), sequenceInfo.orientedReadIds.size());
    }
    for(const SequenceInfo& sequenceInfo: leftFixedSequencesTable) {
        leftFixedSequences.emplace_back(*(sequenceInfo.sequencePointer), sequenceInfo.orientedReadIds.size());
    }
    for(const SequenceInfo& sequenceInfo: rightFixedSequencesTable) {
        rightFixedSequences.emplace_back(*(sequenceInfo.sequencePointer), sequenceInfo.orientedReadIds.size());
    }

    if(html) {
        theseusWriteFile(fixedSequences, leftFixedSequences, rightFixedSequences,
            "Pericles.fasta");
    }

    vector< vector<AlignedBase> > alignment;
    const bool computeAlignment = bool(html);
    theseus(fixedSequences, leftFixedSequences, rightFixedSequences,
        sequence, alignment, computeAlignment);

    if(not html) {
        return;
    }

    SHASTA2_ASSERT(alignment.size() ==
        fixedSequences.size() + leftFixedSequences.size() + rightFixedSequences.size());


    html << "<h4>Alignment</h4>"
        "<table><tr>"
        "<th>Sequence<br>id"
        "<th class=left>Alignment";

    uint64_t sequenceId = 0;
    for(const vector<AlignedBase>& alignmentRow: alignment) {
        html <<
            "<tr>"
            "<td class=centered>" << sequenceId++ <<
            "<td style='font-family:monospace'>";
        std::ranges::copy(alignmentRow, ostream_iterator<AlignedBase>(html));
    }
    html << "</table>";


    html << "<h4>Consensus</h4>"
        "<div style='font-family:monospace'>";
    std::ranges::copy(sequence, ostream_iterator<Base>(html));
    html << "</div>";
}
