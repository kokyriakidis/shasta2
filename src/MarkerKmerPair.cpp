#include "MarkerKmerPair.hpp"
#include "abpoaWrapper.hpp"
#include "MarkerKmers.hpp"
#include "orderPairs.hpp"
#include "Reads.hpp"
using namespace shasta;



MarkerKmerPair::MarkerKmerPair(
    const MarkerKmers& markerKmers,
    const Kmer& kmer0,
    const Kmer& kmer1) :
    kmer0(kmer0),
    kmer1(kmer1)
{
    getMarkerInfos(markerKmers);
    gatherCommonOrientedReads(markerKmers.markers);
    gatherSequences(markerKmers.reads);
    rankSequences();
}



void MarkerKmerPair::getMarkerInfos(const MarkerKmers& markerKmers)
{
    markerKmers.get(kmer0, markerInfos0);
    markerKmers.get(kmer1, markerInfos1);
}



void MarkerKmerPair::gatherCommonOrientedReads(const Markers& markers)
{
    const uint32_t kHalf = uint32_t(markers.k / 2);

    auto it0 = markerInfos0.begin();
    const auto end0 = markerInfos0.end();

    auto it1 = markerInfos1.begin();
    const auto end1 = markerInfos1.end();

    while((it0 != end0) and (it1!=end1)) {

        if(it0->orientedReadId < it1->orientedReadId) {
            ++it0;
        }

        else if(it1->orientedReadId < it0->orientedReadId) {
            ++it1;
        }

        else {

            // We found a common oriented read.
            // If the ordinals are in the correct order, store it.

            if(it0->ordinal < it1->ordinal) {
                const OrientedReadId orientedReadId = it0->orientedReadId;
                const auto orientedReadMarkers = markers[orientedReadId.getValue()];

                commonOrientedReads.emplace_back();
                CommonOrientedRead& commonOrientedRead = commonOrientedReads.back();
                commonOrientedRead.orientedReadId = orientedReadId;
                commonOrientedRead.ordinal0 = it0->ordinal;
                commonOrientedRead.ordinal1 = it1->ordinal;
                commonOrientedRead.position0 = orientedReadMarkers[it0->ordinal].position + kHalf;
                commonOrientedRead.position1 = orientedReadMarkers[it1->ordinal].position + kHalf;
            }

            ++it0;
            ++it1;
        }
    }
}



void MarkerKmerPair::gatherSequences(const Reads& reads)
{
    vector<Base> sequence;

    // Loop over all reads.
    for(uint64_t i=0; i<commonOrientedReads.size(); i++) {
        CommonOrientedRead& commonOrientedRead = commonOrientedReads[i];

        // Get its sequence.
        commonOrientedRead.getSequence(reads, sequence);

        // Find it in the sequence map, adding it if necessary.
        auto it = sequenceMap.find(sequence);
        if(it == sequenceMap.end()) {
            tie(it, ignore) = sequenceMap.insert(make_pair(sequence, SequenceInfo()));
        }
        it->second.orientedReadIndexes.push_back(i);
        commonOrientedRead.sequenceMapIterator = it;
    }
}



void MarkerKmerPair::CommonOrientedRead::getSequence(
    const Reads& reads,
    vector<Base>& sequence) const
{
    sequence.clear();
    for(uint32_t position=position0; position!=position1; position++) {
        sequence.push_back(reads.getOrientedReadBase(orientedReadId, position));
    }
}



void MarkerKmerPair::rankSequences()
{
    // Gather iterators pointing to the sequences and their coverage.
    vector< pair<SequenceMap::iterator, uint64_t> > sequenceTable;
    for(auto it=sequenceMap.begin(); it!=sequenceMap.end(); ++it) {
        sequenceTable.push_back(make_pair(it, it->second.coverage()));
    }

    // Sort them by decreasing coverage.
    std::ranges::sort(sequenceTable, OrderPairsBySecondOnlyGreater<SequenceMap::iterator, uint64_t>());

    // Store the ranks.
    for(uint64_t rank=0; rank<sequenceTable.size(); rank++) {
        const auto it = sequenceTable[rank].first;
        it->second.rank = rank;
        sequencesByRank.push_back(it);
    }
}



void MarkerKmerPair::writeSummary(ostream& html, uint64_t k) const
{
    html <<
        "<h3>Marker k-mer pair summary</h3>"
        "<table>";

    html <<
        "<tr><th class=left>Left k-mer<td class=centered style='font-family:monospace'>";
    kmer0.write(html, k);

    html <<
        "<tr><th class=left>Right k-mer<td class=centered style='font-family:monospace'>";
    kmer1.write(html, k);

    html <<
        "<tr><th class=left>Left k-mer coverage<td class=centered>" << markerInfos0.size() <<
        "<tr><th class=left>Right k-mer coverage<td class=centered>" << markerInfos1.size() <<
        "<tr><th class=left>Common coverage<td class=centered>" << commonOrientedReads.size() <<
        "</table>";

}



void MarkerKmerPair::writeSequences(ostream& html) const
{
    html <<
        "<h3>Distinct oriented read sequences of this marker pair</h3>"
        "<table><tr>"
        "<tr><th>Rank<th>Coverage<th>Length<th class=left>Sequence";

   for(const auto& it: sequencesByRank) {
       const vector<Base>& sequence = it->first;
       const SequenceInfo& sequenceInfo = it->second;
       html <<
           "<tr>"
           "<td class=centered>" << sequenceInfo.rank <<
           "<td class=centered>" << sequenceInfo.coverage() <<
           "<td class=centered>" << sequence.size() <<
           "<td class=left style='font-family:monospace'>";
       std::ranges::copy(sequence, ostream_iterator<Base>(html));

   }

    html << "</table>";
}



void MarkerKmerPair::writeCommonOrientedReads(ostream& html) const
{
    html <<
        "<h3>Common oriented reads</h3>"
        "<table><tr>"
        "<th>Oriented<br>read id"
        "<th>Left<br>ordinal"
        "<th>Right<br>ordinal"
        "<th>Ordinal<br>offset"
        "<th>Left<br>position"
        "<th>Right<br>position"
        "<th>Position<br>offset"
        "<th>Sequence<br>rank"
        "<th>Sequence<br>coverage"
        "<th>Sequence<br>length";
    for(const CommonOrientedRead& commonOrientedRead: commonOrientedReads) {
        const auto& p = *(commonOrientedRead.sequenceMapIterator);
        const vector<Base>& sequence = p.first;
        const SequenceInfo& sequenceInfo = p.second;

        html <<
            "<tr>"
            "<td class=centered>" << commonOrientedRead.orientedReadId <<
            "<td class=centered>" << commonOrientedRead.ordinal0 <<
            "<td class=centered>" << commonOrientedRead.ordinal1 <<
            "<td class=centered>" << commonOrientedRead.ordinalOffset() <<
            "<td class=centered>" << commonOrientedRead.position0 <<
            "<td class=centered>" << commonOrientedRead.position1 <<
            "<td class=centered>" << commonOrientedRead.positionOffset() <<
            "<td class=centered>" << sequenceInfo.rank <<
            "<td class=centered>" << sequenceInfo.coverage() <<
            "<td class=centered>" << sequence.size();
    }
    html << "</table>";

}



void MarkerKmerPair::writeAlignment(ostream& html) const
{
    // Use abpoa to align the sequences.
    vector< vector<Base> > sequences;
    for(const CommonOrientedRead& commonOrientedRead: commonOrientedReads) {
        sequences.push_back(commonOrientedRead.sequenceMapIterator->first);
    }

    // Compute the alignment.
    vector< pair<Base, uint64_t> > consensus;
    vector< vector<AlignedBase> > alignment;
    vector<AlignedBase> alignedConsensus;
    abpoa(sequences, consensus, alignment, alignedConsensus, true);

    // Write the alignment.
    html <<
        "<h3>Alignment</h3>"
        "<table>"
        "<tr><th class=left>OrientedReadId"
        "<th class=left>Sequence<br>length"
        "<th class=left>Aligned sequence";

    for(uint64_t i=0; i<alignment.size(); i++) {
        const vector<AlignedBase>& alignmentRow = alignment[i];

        html << "<tr><th>" << commonOrientedReads[i].orientedReadId <<
            "<td class=centered>" << sequences[i].size() <<
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

    html << "<tr><th>Consensus<td class=centered>" << consensus.size() <<
        "<td style='font-family:monospace;background-color:LightCyan;white-space:nowrap'>";

    uint64_t position = 0;
    for(uint64_t i=0; i<alignedConsensus.size(); i++) {
        const AlignedBase b = alignedConsensus[i];

        if(not b.isGap()) {
            const uint64_t coverage = consensus[position].second;
            html << "<span title='Consensus position " << position << " coverage " << coverage << "'>";
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
            char c;
            if(coverage < 10) {
                c = char(coverage + '0');
            } else if(coverage < 36) {
                c = char(coverage - 10 + 'A');
            } else {
                c = '*';
            }

            html << "<span title='Consensus position " << position << " coverage " << coverage << "'";
            if(coverage < commonOrientedReads.size()) {
                html << " style='background-color:Pink'";
            }
            html << ">";
            html << c << "</span>";

            ++position;
        }
    }

    html << "</table>";



    // Write the consensus.
    html <<
        "<h3>Consensus</h3>"
        "<table>"
        "<tr><th class=left>Consensus sequence length<td class=left>" << consensus.size() <<
        "<tr><th class=left>Consensus sequence"
        "<td style='font-family:monospace'>";

    for(uint64_t position=0; position<consensus.size(); position++) {
        const auto& p = consensus[position];
        const Base b = p.first;
        const uint64_t coverage = p.second;
        html << "<span title='Consensus position " << position << " coverage " << coverage << "'";
        if(coverage < commonOrientedReads.size()) {
            html << "style='background-color:Pink'";
        }
        html << ">" << b << "</span>";
    }

    html <<
        "<tr><th class=left >Coverage"
        "<td style='font-family:monospace'>";

    std::map<char, uint64_t> coverageLegend;

    bool hasHighCoverage = false;
    for(uint64_t position=0; position<consensus.size(); position++) {
        const auto& p = consensus[position];
        const uint64_t coverage = p.second;

        char c;
        if(coverage < 10) {
            c = char(coverage + '0');
        } else if(coverage < 36) {
            c = char(coverage - 10 + 'A');
        } else {
            c = '*';
            hasHighCoverage = true;
        }
        if(c != '*') {
            coverageLegend.insert(make_pair(c, coverage));
        }

        if(coverage < commonOrientedReads.size()) {
            html << "<span style='background-color:Pink' "
                "title='Consensus position " << position << " coverage " << coverage << "'>";
        }

        html << c;

        if(coverage < commonOrientedReads.size()) {
            html << "</span>";
        }
    }

    html << "</table>";



    // Write the coverage legend.
    html << "<p><table><tr><th>Symbol<th>Coverage";
    for(const auto& p: coverageLegend) {
        html << "<tr><td class=centered>" << p.first << "<td class=centered>" << p.second;
    }
    if(hasHighCoverage) {
        html << "<tr><td class=centered>*<td class=centered>&ge;36";
    }
    html << "</table>";
}



void MarkerKmerPair::writePairAlignmentDistances(ostream& html) const
{
    const uint64_t n = sequencesByRank.size();
    vector< vector<uint64_t> > editDistances(n, vector<uint64_t>(n, 0));

    vector< vector<Base> > sequences(2);
    vector< pair<Base, uint64_t> > consensus;
    vector< vector<AlignedBase> > alignment;
    vector<AlignedBase> alignedConsensus;

    // Loop over pairs of sequences.
    for(uint64_t i1=1; i1<n; i1++) {
        const vector<Base>& sequence1 = sequencesByRank[i1]->first;
        sequences[1] = sequence1;

        for(uint64_t i0=0; i0<i1; i0++) {
            const vector<Base>& sequence0 = sequencesByRank[i0]->first;
            sequences[0] = sequence0;

            // Align them.
            abpoa(sequences, consensus, alignment, alignedConsensus, true);

            // Compute edit distances.
            uint64_t editDistance = 0;
            for(uint64_t position=0; position<alignment[0].size(); position++) {
                if(alignment[0][position] != alignment[1][position]) {
                    ++editDistance;
                }
            }
            editDistances[i0][i1] = editDistance;
            editDistances[i1][i0] = editDistance;
        }
    }


    // Write out the edit distances.
    html <<
        "<h3>Edit distances between the sequences</h3>"
        "<p><table><tr><th>Ranks";
    for(uint64_t i=0; i<n; i++) {
        html << "<th>" << i;
    }

    for(uint64_t i0=0; i0<n; i0++) {
        html << "<tr><th>" << i0;
        for(uint64_t i1=0; i1<n; i1++) {
            html << "<td class=centered>" << editDistances[i0][i1];
        }
    }

    html << "</table>";
}
