// Shasta.
#include "LocalAssembly6.hpp"
#include "Anchor.hpp"
#include "Base.hpp"
#include "orderPairs.hpp"
#include "Reads.hpp"
#include "theseusWrapper.hpp"
using namespace shasta2;

// Standard library.
#include <cmath>
#include <map>
#include <set>



LocalAssembly6::LocalAssembly6(
    const Anchors& anchors,
    AnchorId anchorIdA,
    AnchorId anchorIdB,
    ostream& html,
    const vector<OrientedReadId>& orientedReadIds) :
    anchors(anchors),
    anchorIdA(anchorIdA),
    anchorIdB(anchorIdB),
    html(html)
{
    if(html) {
        html << "<h3>Local assembly between " <<
            anchorIdToString(anchorIdA) << " and " <<
            anchorIdToString(anchorIdB) << "</h3>";
    }

    gatherOrientedReads(orientedReadIds);
    removeOutliers();
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

        // If on both anchors and negative offset, this OrientedReadId cannot be used.
        if(isOnA and isOnB and (itA->position > itB->position)) {
            continue;
        }

        // We can use this OrientedReadId for assembly.

        // Create an OrientedReadInfo.
        OrientedReadInfo& info = orientedReadInfos.emplace_back();
        info.orientedReadId = orientedReadId;
        if(isOnA) {
            info.positionA = itA->position;
        }
        if(isOnB) {
            info.positionB = itB->position;
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
            "<td class=left style='font-family:monospace;white-space: nowrap'>";
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
    vector<AlignedBase> alignedConsensus;
    vector<uint64_t> coverage;
    theseus(fixedSequences, leftFixedSequences, rightFixedSequences,
        sequence, alignedConsensus, coverage, alignment, computeAlignment);

    if(not html) {
        return;
    }

    SHASTA2_ASSERT(alignment.size() ==
        fixedSequences.size() + leftFixedSequences.size() + rightFixedSequences.size());



    // Write the alignment.
    html <<
        "<h4>Alignment</h4>"
        "<table><tr>"
        "<th>Sequence<br>id"
        "<th>Fixed<br>on A"
        "<th>Fixed<br>on B"
        "<th>Length"
        "<th>Coverage"
        "<th class=left>Alignment";

    for(uint64_t sequenceId=0; sequenceId<alignment.size(); sequenceId++) {
        const vector<AlignedBase>& alignmentRow = alignment[sequenceId];
        bool fixedOnA = false;
        bool fixedOnB = false;
        const SequenceInfo& sequenceInfo = getSequenceInfo(sequenceId, fixedOnA, fixedOnB);
        html <<
            "<tr>"
            "<td class=centered>" << sequenceInfo.id <<
            "<td class=centered>" << (fixedOnA ? "&check;" : "") <<
            "<td class=centered>" << (fixedOnB ? "&check;" : "") <<
            "<td class=centered>" << sequenceInfo.sequencePointer->size() <<
            "<td class=centered>" << sequenceInfo.orientedReadIds.size() <<
            "<td class=left style='font-family:monospace;white-space:nowrap'>";
        for(uint64_t position=0; position<alignmentRow.size(); position++) {
            const AlignedBase base = alignmentRow[position];
            const bool isMismatch = (base != alignedConsensus[position]);
            if(isMismatch) {
                html << "<span style='background-color:Pink'>";
            }
            html << base;
            if(isMismatch) {
                html << "</span>";
            }
        }
    }



    // Write a line with aligned consensus.
    html <<
        "<tr>"
        "<th colspan=3 class=left>Aligned consensus"
        "<td class=centered>" << sequence.size() <<
        "<td>"
        "<td class=left style='font-family:monospace;white-space:nowrap'>";
    uint64_t nonAlignedPosition = 0;
    for(uint64_t position=0; position<alignedConsensus.size(); position++) {
        const AlignedBase base = alignedConsensus[position];
        if(base.isGap()) {
            html << "-";
        } else {
            const uint64_t c = coverage[nonAlignedPosition];
            html << "<span title='";
            html << "Position " << nonAlignedPosition << " coverage " << c;
            html << "'";
            html << ">";
            html << base;
            html << "</span>";
             ++nonAlignedPosition;
        }
    }
    SHASTA2_ASSERT(nonAlignedPosition == sequence.size());



    // Write a line with aligned consensus coverage.
    html <<
        "<tr>"
        "<th colspan=3 class=left>Coverage"
        "<td>"
        "<td>"
        "<td class=left style='font-family:monospace;white-space:nowrap'>";
    nonAlignedPosition = 0;
    std::map<char, uint64_t> coverageLegend;
    for(uint64_t position=0; position<alignedConsensus.size(); position++) {
        if(alignedConsensus[position].isGap()) {
            html << "-";
        } else {
            const uint64_t c = coverage[nonAlignedPosition];
            char coverageCharacter = ' ';
            if(c < 10) {
                coverageCharacter = char(c - '0');
            } else if(c < 36) {
                coverageCharacter = char(c - 10 + 'A');
            } else {
                coverageCharacter = '*';
            }
            html << coverageCharacter;
            coverageLegend[coverageCharacter] = c;
            ++nonAlignedPosition;
        }
    }
    SHASTA2_ASSERT(nonAlignedPosition == sequence.size());

    html << "</table>";



    // Write the consensus.
    html <<
        "<h4>Consensus</h4>"
        "<table>"
        "<tr><th class=left>Length<td class=left>" << sequence.size() <<
        "<tr><th class=left>Sequence<td class=left style='font-family:monospace;white-space:nowrap'>";
    for(uint64_t position=0; position<sequence.size(); position++) {
        html << "<span title='Position " << position <<
            " coverage " << coverage[position] << "'";
        html << ">";
        html << sequence[position];
        html << "</span>";

    }
    html <<
        "<tr><th class=left>Coverage<td class=left style='font-family:monospace;white-space:nowrap'>";
    for(uint64_t position=0; position<sequence.size(); position++) {
        const uint64_t c = coverage[position];
        char coverageCharacter = ' ';
        if(c < 10) {
            coverageCharacter = char(c - '0');
        } else if(c < 36) {
            coverageCharacter = char(c - 10 + 'A');
        } else {
            coverageCharacter = '*';
        }
        html << coverageCharacter;
    }
    html << "</table>";



    // Write coverage legend.
    html << "<h4>Coverage legend</h4>"
        "<table><tr><td class=centered>Character<td class=centered>Coverage";
    for(const auto& [character, coverage]: coverageLegend) {
        html << "<tr><td class=centered>" << character <<
            "<td class=centered>" << coverage;
    }
    html << "</table>";



}



// If the offsets of oriented reads constrained at both A and B
// are too different, remove the outliers.
void LocalAssembly6::removeOutliers()
{
    // Gather the offsets.
    vector<pair<uint64_t, uint64_t > > offsetTable; // (index in orientedReadInfos, offset).
    for(uint64_t i=0; i<orientedReadInfos.size(); i++) {
        const OrientedReadInfo& orientedReadInfo = orientedReadInfos[i];
        if(orientedReadInfo.isOnBothAnchors()) {
            const uint64_t offset = orientedReadInfo.positionOffsetAB();
            offsetTable.push_back({i, offset});
        }
    }
    sort(offsetTable.begin(), offsetTable.end(), OrderPairsBySecondOnly<uint64_t, uint64_t>());

    // Find places where there is an unreasonably jump in the offset.
    vector<uint64_t> violations(1, 0);
    for(uint64_t i1=1; i1<offsetTable.size(); i1++) {
        const uint64_t i0 = i1 - 1;
        const uint64_t offset0 = offsetTable[i0].second;
        const uint64_t offset1 = offsetTable[i1].second;
        if(not checkOffsets(offset0, offset1)) {
            violations.push_back(i1);
            // cout << "Violation " << offset0 << " " << offset1 << " " << i1 << endl;
        }
    }
    violations.push_back(offsetTable.size());

    // If no violations were found, keep all the OrientedReadInfos.
    // This is the most common case.
    if(violations.size() == 2) {
        return;
    }

#if 0
    cout << "violations vector ";
    std::ranges::copy(violations, ostream_iterator<uint64_t>(cout, " "));
    cout << endl;
#endif

    // Find the largest interval between violations.
    uint64_t keepBegin = 0;
    uint64_t keepEnd = 0;
    for(uint64_t i1=1; i1<violations.size(); i1++) {
        const uint64_t i0 = i1 - 1;
        const uint64_t violation0 = violations[i0];
        const uint64_t violation1 = violations[i1];
        if(violation1 - violation0 > keepEnd - keepBegin) {
            keepBegin = violation0;
            keepEnd = violation1;
        }
    }
    // cout << "keepBegin " << keepBegin << ", keepEnd " << keepEnd << endl;

    // Only keep OrientedReadInfos that are at positions [keepBegin, keepEnd)
    // in the offset table.
    std::set<uint64_t> discard;
    for(uint64_t i=0; i<keepBegin; i++) {
        const uint64_t j = offsetTable[i].first;
        discard.insert(j);
        if(html) {
            html << "<br>Discarding " << orientedReadInfos[j].orientedReadId <<
                " due to inconsistent offsets.";
        }
    }
    for(uint64_t i=keepEnd; i<offsetTable.size(); i++) {
        const uint64_t j = offsetTable[i].first;
        discard.insert(j);
        if(html) {
            html << "<br>Discarding " << orientedReadInfos[j].orientedReadId <<
                " due to inconsistent offsets.";
        }
    }

    vector<OrientedReadInfo> newOrientedReadInfos;
    for(uint64_t i=0; i<orientedReadInfos.size(); i++) {
        if(not discard.contains(i)) {
            newOrientedReadInfos.push_back(orientedReadInfos[i]);
        }
    }
    orientedReadInfos.swap(newOrientedReadInfos);
}



bool LocalAssembly6::checkOffsets(uint64_t offset0, uint64_t offset1)
{
    // EXPOSE WHEN CODE STABILIZES.
    const double aDrift = 0.02;
    const double bDrift = 100.;

    if(offset1 == offset0) {
        return true;
    }

    SHASTA2_ASSERT(offset1 > offset0);

    const double average = 0.5 * double(offset0 + offset1);
    const uint64_t difference = offset1 - offset0;

    const double acceptableDifference = aDrift * average + bDrift;

    return difference < uint64_t(std::round(acceptableDifference));

}



// The sequences are passed to theseus in the above order.
// This returns the SequenceInfo with a given index in that order.
const LocalAssembly6::SequenceInfo& LocalAssembly6::getSequenceInfo(
    uint64_t sequenceId,
    bool& fixedOnA,
    bool& fixedOnB) const
{
    if(sequenceId < fixedSequencesTable.size()) {
        fixedOnA = true;
        fixedOnB = true;
        return fixedSequencesTable[sequenceId];
    }
    sequenceId -= fixedSequencesTable.size();
    if(sequenceId < leftFixedSequencesTable.size()) {
        fixedOnA = true;
        fixedOnB = false;
        return leftFixedSequencesTable[sequenceId];
    }
    sequenceId -= leftFixedSequencesTable.size();
    if(sequenceId < rightFixedSequencesTable.size()) {
        fixedOnA = false;
        fixedOnB = true;
        return rightFixedSequencesTable[sequenceId];
    } else {
        SHASTA2_ASSERT(0);
    }
}
