// Shasta.
#include "LocalAssembly4.hpp"
#include "Anchor.hpp"
#include "Markers.hpp"
#include "poastaWrapper.hpp"
#include "ReadId.hpp"
#include "Reads.hpp"
using namespace shasta2;

// Standard library.
#include "algorithm.hpp"
#include <map>



LocalAssembly4::LocalAssembly4(
    const Anchors& anchors,
    [[maybe_unused]] uint64_t abpoaMaxLength,
    ostream& html,
    bool /* debug */,
    const AnchorPair& anchorPair,
    const vector<OrientedReadId>& additionalOrientedReadIds) :
    anchors(anchors),
    html(html),
    leftAnchorId(anchorPair.anchorIdA),
    rightAnchorId(anchorPair.anchorIdB)
{
    SHASTA2_ASSERT(std::ranges::is_sorted(additionalOrientedReadIds));
    if(html) {
        writeInput(anchorPair, additionalOrientedReadIds);
    }

    gatherAllOrientedReads(anchorPair, additionalOrientedReadIds);
    if(html) {
        writeAllOrientedReadIds();
    }

    gatherCommonOrientedReads();
    if(html) {
        writeCommonOrientedReads();
    }

    assemble();
    if(html) {
        writeAssembledSequence();
    }
}



void LocalAssembly4::gatherAllOrientedReads(
    const AnchorPair& anchorPair,
    const vector<OrientedReadId>& additionalOrientedReadIds)
{
    std::ranges::set_union(anchorPair.orientedReadIds, additionalOrientedReadIds,
        back_inserter(allOrientedReadIds));

}



// Gather information on oriented reads in the usableOrientedReadIds
// that appear in both the left and right AnchorIds.
void LocalAssembly4::gatherCommonOrientedReads()
{

    const uint32_t kHalf = uint32_t(anchors.k / 2);

    // The common OrientedReadIds must appear in both the
    // left and right anchors (the latter visited after the former).
    // These are the OrientedReadIds in the "complete" AnchorPair.
    const AnchorPair completeAnchorPair(anchors, leftAnchorId, rightAnchorId, false);

    // The common OrientedReadIds are the intersection of the
    // allOrientedReadIds and the OrientedReadIds in the completeAnchorPair.
    vector<OrientedReadId> commonOrientedReadIds;
    std::ranges::set_intersection(allOrientedReadIds, completeAnchorPair.orientedReadIds,
        back_inserter(commonOrientedReadIds));



    // Now we can fill in the commonOrientedReadInfos.
    for(const OrientedReadId orientedReadId: commonOrientedReadIds) {
        const auto orientedReadMarkers = anchors.markers[orientedReadId.getValue()];

        CommonOrientedReadInfo& commonOrientedReadInfo = commonOrientedReadInfos.emplace_back();
        commonOrientedReadInfo.orientedReadId = orientedReadId;

        commonOrientedReadInfo.leftOrdinal = anchors.getOrdinal(leftAnchorId, orientedReadId);
        commonOrientedReadInfo.rightOrdinal = anchors.getOrdinal(rightAnchorId, orientedReadId);

        commonOrientedReadInfo.leftPosition = orientedReadMarkers[commonOrientedReadInfo.leftOrdinal].position + kHalf;
        commonOrientedReadInfo.rightPosition = orientedReadMarkers[commonOrientedReadInfo.rightOrdinal].position + kHalf;
     }
}



void LocalAssembly4::assemble()
{
    SHASTA2_ASSERT(not commonOrientedReadInfos.empty());

    // Gather distinct sequences and their coverage.
    class DistinctSequence {
    public:
        shared_ptr< vector<Base> > sequencePointer;
        const vector<Base>& sequence() const
        {
            return *sequencePointer;
        }
        uint64_t coverage;

        // Order by decreasing coverage.
        bool operator<(const DistinctSequence& that) {
            return coverage > that.coverage;
        }
    };
    vector<DistinctSequence> distinctSequences;
    vector<Base> orientedReadSequence;
    for(const CommonOrientedReadInfo& info: commonOrientedReadInfos) {
        getSequence(info, orientedReadSequence);

        // See if this sequence is already in our DistinctSequences.
        bool done = false;
        for(DistinctSequence& distinctSequence: distinctSequences) {
            if(distinctSequence.sequence() == orientedReadSequence) {
                ++distinctSequence.coverage;
                done = true;
                break;
            }
        }

        // If we did not have it, create it with coverage 1.
        if(not done) {
            DistinctSequence& distinctSequence = distinctSequences.emplace_back();
            distinctSequence.sequencePointer = make_shared< vector<Base> >(orientedReadSequence);
            distinctSequence.coverage = 1;
        }
    }
    sort(distinctSequences.begin(), distinctSequences.end());



    // Write the DistinctSequences.
    if(html) {
        html <<
            "<h3>Distinct sequences</h3>"
            "<table><tr>"
            "<th>Coverage<th>Length<th>Sequence";
        for(const DistinctSequence& distinctSequence: distinctSequences) {
            const vector<Base>& sequence = distinctSequence.sequence();
            html <<
                "<tr>"
                "<td class=centered>" << distinctSequence.coverage <<
                "<td class=centered>" << sequence.size() <<
                "<td class=left style='font-family:monospace;white-space:nowrap''>";
            std::ranges::copy(sequence, ostream_iterator<Base>(html));
        }
        html << "</table>";
    }



    // If there is a dominant sequence, use it as the consensus.
    const uint64_t maximumCoverage = distinctSequences.front().coverage;
    const uint64_t totalCoverage = commonOrientedReadInfos.size();
    if(maximumCoverage > totalCoverage/2) {
        sequence = distinctSequences.front().sequence();
        coverage.resize(sequence.size(), maximumCoverage);
        if(html) {
            html << "<br>The dominant sequence with coverage " << maximumCoverage <<
                " was used as the consensus.";
            if(maximumCoverage < totalCoverage) {
                html << "<br>Store base coverages are lower bounds.";
            }
        }
        return;
    }



    // If getting here, we have to run the multiple sequence alignment of the distinct sequences.
    vector< pair<vector<Base>, uint64_t> > sequencesWithCoverage;
    for(const DistinctSequence& distinctSequence: distinctSequences) {
         sequencesWithCoverage.emplace_back(distinctSequence.sequence(), distinctSequence.coverage);
    }
    vector< pair<Base, uint64_t> > consensus;
    vector< vector<AlignedBase> > alignment;
    vector<AlignedBase> alignedConsensus;
    poasta(sequencesWithCoverage, consensus, alignment, alignedConsensus);

    // Store assembled sequence
    for(const auto& [base, baseCoverage]: consensus) {
        sequence.push_back(base);
        coverage.push_back(baseCoverage);
    }



    // Write the multiple sequence alignment of the distinct sequences.
    if(html) {
        html <<
            "<h3>Multiple sequence alignment of the distinct sequences</h3>"
            "<table><tr>"
            "<th>Coverage<th>Length<th>Aligned sequence";
        for(uint64_t i=0; i<distinctSequences.size(); i++) {
            const DistinctSequence& distinctSequence = distinctSequences[i];
            const vector<Base>& sequence = distinctSequence.sequence();
            const vector<AlignedBase>& alignmentRow = alignment[i];
            html <<
                "<tr>"
                "<td class=centered>" << distinctSequence.coverage <<
                "<td class=centered>" << sequence.size() <<
                "<td class=left style='font-family:monospace;white-space:nowrap''>";
            std::ranges::copy(alignmentRow, ostream_iterator<AlignedBase>(html));
        }
        html << "</table>";
    }
}



// Get the sequence a CommonOrientedReadInfo contributes to the assembly.
void LocalAssembly4::getSequence(const CommonOrientedReadInfo& info, vector<Base>& sequence) const
{
    sequence.clear();
    for(uint32_t position=info.leftPosition; position<info.rightPosition; position++) {
        sequence.push_back(anchors.reads.getOrientedReadBase(info.orientedReadId, position));
    }
}



void LocalAssembly4::writeInput(
    const AnchorPair& anchorPair,
    const vector<OrientedReadId>& additionalOrientedReadIds) const
{
    html <<
        "<table>"
        "<tr><th class=left>Left anchor<td class=centered>" << anchorIdToString(leftAnchorId) <<
        "<tr><th class=left>Right anchor<td class=centered>" << anchorIdToString(rightAnchorId) <<
        "</table>";
    html << "<h3>OrientedReadIds in the AnchorPair</h3><table>";
    for(const OrientedReadId orientedReadId: anchorPair.orientedReadIds) {
        html << "<tr><td class=centered>" << orientedReadId;
    }
    html << "</table>" << anchorPair.orientedReadIds.size() << " oriented reads in the AnchorPair.";

    html << "<h3>Additional OrientedReadIds</h3>"
        "These are additional OrientedReadIds that can be used in this assembly."
        "<table>";
    for(const OrientedReadId orientedReadId: additionalOrientedReadIds) {
        html << "<tr><td class=centered>" << orientedReadId;
    }
    html << "</table>" << additionalOrientedReadIds.size() << " additional oriented reads.";

}


void LocalAssembly4::writeAllOrientedReadIds() const
{
    html << "<h3>All OrientedReadIds</h3>"
        "These are the union of the OrientedReadIds in the AnchorPair "
        "and the additional OrientedReadIds"
        "<table>";
    for(const OrientedReadId orientedReadId: allOrientedReadIds) {
        html << "<tr><td class=centered>" << orientedReadId;
    }
    html << "</table>" << allOrientedReadIds.size() << " oriented reads.";
}



void LocalAssembly4::writeCommonOrientedReads() const
{
    html <<
        "<h3>Common oriented reads</h3>"
        "These are OrientedReadIds from the above list that are present in both the "
        "left and right anchors, and that visit the right anchor "
        "after (not necessarily immediately after) they visit the left anchor."
        "<table>"
        "<tr>"
        "<th>OrientedReadId"
        "<th>Left<br>ordinal<th>Right<br>ordinal<th>Ordinal<br>offset"
        "<th>Left<br>position<th>Right<br>position<th>Position<br>offset"
        "<th>Sequence";

    vector<Base> sequence;
    for(const CommonOrientedReadInfo& info: commonOrientedReadInfos) {
        getSequence(info, sequence);
        html <<
            "<tr>"
            "<td class=centered>" << info.orientedReadId <<
            "<td class=centered>" << info.leftOrdinal <<
            "<td class=centered>" << info.rightOrdinal <<
            "<td class=centered>" << info.ordinalOffset() <<
            "<td class=centered>" << info.leftPosition <<
            "<td class=centered>" << info.rightPosition <<
            "<td class=centered>" << info.positionOffset() <<
            "<td class=left style='font-family:monospace;white-space: nowrap''>";
        std::ranges::copy(sequence, ostream_iterator<Base>(html));
    }
    html << "</table>" << commonOrientedReadInfos.size() << " common oriented reads.";

}



void LocalAssembly4::writeAssembledSequence() const
{
    html <<
        "<h3>Consensus</h3>"
        "<p><span style='font-family:monospace;white-space:nowrap''>"
        ">LocalAssembly " << sequence.size() <<
        "<br>";
    std::ranges::copy(sequence, ostream_iterator<Base>(html));
    html << "</span><br><br>";

    html <<
        "<table>"
        "<tr><th class=left>Consensus sequence length<td class=left>" << sequence.size() <<
        "<tr><th class=left>Consensus sequence"
        "<td style='font-family:monospace;white-space:nowrap''>";

    for(uint64_t position=0; position<sequence.size(); position++) {
        const Base b = sequence[position];
        html << "<span title='" << position << "'>" << b << "</span>";
    }

    html <<
        "<tr><th class=left >Coverage"
        "<td style='font-family:monospace;white-space:nowrap''>";

    std::map<char, uint64_t> coverageLegend;

    for(uint64_t position=0; position<sequence.size(); position++) {
        const uint64_t coverageThisPosition = coverage[position];
        const char c = (coverageThisPosition < 10) ? char(coverageThisPosition + '0') : char(coverageThisPosition - 10 + 'A');
        coverageLegend.insert(make_pair(c, coverageThisPosition));

        if(coverageThisPosition < commonOrientedReadInfos.size()) {
            html << "<span style='background-color:Pink'>";
        }

        html << c;

        if(coverageThisPosition < commonOrientedReadInfos.size()) {
            html << "</span>";
        }
    }

    html << "</table><br><br>";

    // Write the coverage legend.
    html << "<p><table><tr><th>Symbol<th>Coverage";
    for(const auto& p: coverageLegend) {
        html << "<tr><td class=centered>" << p.first << "<td class=centered>" << p.second;
    }
    html << "</table>";


}
