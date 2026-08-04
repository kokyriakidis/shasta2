#include "theseusWrapper.hpp"
#include "Base.hpp"
#include "msa1.hpp"
#include "SHASTA2_ASSERT.hpp"
using namespace shasta2;

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-conversion"
#include "theseus/theseus_msa_aligner.h"
#pragma GCC diagnostic pop

#include <boost/tokenizer.hpp>

#include "fstream.hpp"
#include "iterator.hpp"
#include <sstream>



namespace shasta2 {

    // The penalties and heuristics used for all theseus alignments here.
    // Pericles defaults.
    static theseus::Penalties theseusPenalties()
    {
        const int match = 0;
        const int mismatch = 2;
        const int gapo = 3;
        const int gape = 1;
        return theseus::Penalties(match, mismatch, gapo, gape);
    }



    // Feed a group of sequences, already converted to strings, to a theseus
    // aligner, using the flags appropriate for each group.
    // The first sequence fixed on both sides must already have been used to
    // construct the aligner, so it is skipped here.
    static void theseusAlignGroups(
        theseus::TheseusMSA& aligner,
        const vector< pair<string, uint64_t> >& fixedSequences,
        const vector< pair<string, uint64_t> >& leftFixedSequences,
        const vector< pair<string, uint64_t> >& rightFixedSequences)
    {
        // The remaining sequences fixed on both sides.
        for(uint64_t i=1; i<fixedSequences.size(); i++) {
            const auto& [sequence, weight] = fixedSequences[i];
            const bool densityDrop = false;
            const bool lagPruning = false;
            aligner.align(sequence, int(weight), false, false, densityDrop, lagPruning);
        }

        // The sequences fixed on the left only.
        for(const auto& [sequence, weight]: leftFixedSequences) {
            const bool densityDrop = true;
            const bool lagPruning = false;
            aligner.align(sequence, int(weight), false, true, densityDrop, lagPruning);
        }

        // The sequences fixed on the right only.
        for(const auto& [sequence, weight]: rightFixedSequences) {
            const bool densityDrop = true;
            const bool lagPruning = false;
            aligner.align(sequence, int(weight), true, true, densityDrop, lagPruning);
        }
    }



    // Extract the alignment rows from a theseus aligner, as strings.
    // print_as_msa writes a header line for each sequence, then its alignment
    // row, and finally a row for the aligned consensus, which is dropped here.
    static void theseusAlignmentRows(
        theseus::TheseusMSA& aligner,
        uint64_t expectedRowCount,
        vector<string>& rows)
    {
        rows.clear();
        std::ostringstream s;
        aligner.print_as_msa(s);
        const string alignmentString = s.str();

        boost::tokenizer< boost::char_separator<char> > tokenizer(
            alignmentString, boost::char_separator<char>("\n"));
        for(const string& line: tokenizer) {
            SHASTA2_ASSERT(not line.empty());
            if(line[0] != '>') {
                rows.push_back(line);
            }
        }

        // Drop the aligned consensus row.
        SHASTA2_ASSERT(not rows.empty());
        rows.pop_back();
        SHASTA2_ASSERT(rows.size() == expectedRowCount);
    }
}



// See theseusWrapper.hpp for comments.
void shasta2::theseusExtended(
    const vector< pair<vector<ExtendedBase>, uint64_t> >& fixedSequences,
    const vector< pair<vector<ExtendedBase>, uint64_t> >& leftFixedSequences,
    const vector< pair<vector<ExtendedBase>, uint64_t> >& rightFixedSequences,
    vector< vector<ExtendedBase> >& alignment)
{
    // Convert to the ACGTacgt string form theseus works on.
    // Lower case marks a poly symbol, and theseus compares characters directly,
    // so the distinction survives the round trip.
    const auto toStrings = [](
        const vector< pair<vector<ExtendedBase>, uint64_t> >& in,
        vector< pair<string, uint64_t> >& out)
    {
        out.clear();
        for(const auto& [sequence, weight]: in) {
            out.push_back(make_pair(toString(sequence), weight));
        }
    };
    vector< pair<string, uint64_t> > fixedStrings;
    vector< pair<string, uint64_t> > leftFixedStrings;
    vector< pair<string, uint64_t> > rightFixedStrings;
    toStrings(fixedSequences, fixedStrings);
    toStrings(leftFixedSequences, leftFixedStrings);
    toStrings(rightFixedSequences, rightFixedStrings);

    SHASTA2_ASSERT(not fixedStrings.empty());
    const theseus::Penalties penalties = theseusPenalties();
    theseus::Heuristics heuristics;
    theseus::TheseusMSA aligner(penalties, heuristics,
        fixedStrings.front().first, int(fixedStrings.front().second), false);

    theseusAlignGroups(aligner, fixedStrings, leftFixedStrings, rightFixedStrings);

    vector<string> rows;
    theseusAlignmentRows(aligner,
        fixedSequences.size() + leftFixedSequences.size() + rightFixedSequences.size(),
        rows);

    // Parse back with ExtendedBase, NOT AlignedBase. AlignedBase::fromCharacter
    // maps 'a' to A, which would silently discard the poly distinction.
    alignment.clear();
    for(const string& row: rows) {
        alignment.push_back(vectorOfExtendedBasesFromString(row));
    }

    // All rows must have the same length.
    for(const vector<ExtendedBase>& row: alignment) {
        SHASTA2_ASSERT(row.size() == alignment.front().size());
    }
}



void shasta2::theseus(

    // The input sequences fixed on both sides, with their coverage.
    // They are passed to theseus in this order.
    const vector< pair<vector<Base>, uint64_t> >& fixedSequences,

    // The input sequences fixed on the left only, with their coverage.
    const vector< pair<vector<Base>, uint64_t> >& leftFixedSequences,

    // The input sequences fixed on the right only, with their coverage.
    const vector< pair<vector<Base>, uint64_t> >& rightFixedSequences,

    // The consensus sequence and its coverage.
    vector< pair<Base, uint64_t> >& consensus,

    // The computed alignment.
    // Each element of the vector correspond to one of the input sequences,
    // in the same order.
    // These all have the same length, which equals the length of the aligned consensus.
    vector< vector<AlignedBase> >& alignment,

    // The aligned consensus.
    vector<AlignedBase>& alignedConsensus,

    // Consensus and alignedConsensus are always computed.
    // Alignment is only computed if this set to true.
    bool computeAlignment
)
{


    // Convert to strings.
    const auto toStrings = [](
        const vector< pair<vector<Base>, uint64_t> >& in,
        vector< pair<string, uint64_t> >& out)
    {
        out.clear();
        for(const auto& [sequence, weight]: in) {
            out.push_back(make_pair(toString(sequence), weight));
        }
    };
    vector< pair<string, uint64_t> > fixedStrings;
    vector< pair<string, uint64_t> > leftFixedStrings;
    vector< pair<string, uint64_t> > rightFixedStrings;
    toStrings(fixedSequences, fixedStrings);
    toStrings(leftFixedSequences, leftFixedStrings);
    toStrings(rightFixedSequences, rightFixedStrings);

    // Create the theseus aligner, passing in the first sequence fixed on both sides.
    SHASTA2_ASSERT(not fixedStrings.empty());
    const theseus::Penalties penalties = theseusPenalties();
    theseus::Heuristics heuristics;
    theseus::TheseusMSA aligner(penalties, heuristics,
        fixedStrings.front().first, int(fixedStrings.front().second), false);

    // Pass the remaining sequences to the aligner.
    theseusAlignGroups(aligner, fixedStrings, leftFixedStrings, rightFixedStrings);

    // Get the consensus.
    vector<int> coverage;
    string consensusString;
    string alignedConsensusString;
    aligner.majority_voting_consensus(coverage, consensusString, alignedConsensusString);
    SHASTA2_ASSERT(coverage.size() == consensusString.size());
    consensus.clear();
    for(uint64_t i=0; i<consensusString.size(); i++) {
        consensus.push_back(make_pair(Base::fromCharacter(consensusString[i]), coverage[i]));
    }
    alignedConsensus = vectorOfAlignedBasesFromString(alignedConsensusString);

    // Also compute the alignment, if requested.
    if(computeAlignment) {
        vector<string> rows;
        theseusAlignmentRows(aligner,
            fixedSequences.size() + leftFixedSequences.size() + rightFixedSequences.size(),
            rows);
        alignment.clear();
        for(const string& row: rows) {
            alignment.push_back(vectorOfAlignedBasesFromString(row));
        }
    }
}



void shasta2::testTheseus()
{
    const vector< pair<vector<Base>, uint64_t> > fixedSequences =
    {
        {vectorOfBasesFromString("TAGGGATTGATAAAAGGCTTTCCAGAAGA"), 5},
        {vectorOfBasesFromString("TAGGGATTCATAAAGGCTTTCCAGAAGA"), 4},
        {vectorOfBasesFromString("TAGGGATTGATAAAAAGGCTTTCCAGAAGA"), 3}
    };

    const vector< pair<vector<Base>, uint64_t> > leftFixedSequences =
    {
        {vectorOfBasesFromString("TAGGGATTGATAAAAGGCTTTCCAGAAGATTTTTTT"), 10}
    };

    const vector< pair<vector<Base>, uint64_t> > rightFixedSequences =
    {
        {vectorOfBasesFromString("GGGGGGTAGGGATTGATAAAAGGCTTTCCAGAAGA"), 10}
    };




    vector< pair<Base, uint64_t> > consensus;
    vector<AlignedBase> alignedConsensus;
    vector< vector<AlignedBase> > alignment;
    theseus(fixedSequences, leftFixedSequences, rightFixedSequences,
        consensus, alignment, alignedConsensus, true);

    cout << "Consensus with coverage:" << endl;
    for(const auto& [base, coverage]: consensus) {
        cout << base;
    }
    cout << endl;



    cout << "Alignment:" << endl;
    for(const vector<AlignedBase>& alignmentRow: alignment) {
        for(const AlignedBase base: alignmentRow) {
            cout << base;
        }
        cout << endl;
    }

    cout << "Aligned consensus:" << endl;
    std::ranges::copy(alignedConsensus, ostream_iterator<AlignedBase>(cout));
    cout << endl;

}



// This writes a file that can be used as input to pericles.
void shasta2::theseusWriteFile(

    // The input sequences fixed on both sides, with their coverage.
    // They are passed to theseus in this order.
    const vector< pair<vector<Base>, uint64_t> >& fixedSequences,

    // The input sequences fixed on the left only, with their coverage.
    const vector< pair<vector<Base>, uint64_t> >& leftFixedSequences,

    // The input sequences fixed on the right only, with their coverage.
    const vector< pair<vector<Base>, uint64_t> >& rightFixedSequences,

    const string& fileName
)
{
    ofstream out(fileName);

    for(const auto& [sequence, weight]: fixedSequences) {
        out << ">0 0 " << weight << " 0 0 \n";
        std::ranges::copy(sequence, ostream_iterator<Base>(out));
        out << "\n";
    }

    for(const auto& [sequence, weight]: leftFixedSequences) {
        out << ">0 1 " << weight << " 1 0\n";
        std::ranges::copy(sequence, ostream_iterator<Base>(out));
        out << "\n";
    }

    for(const auto& [sequence, weight]: rightFixedSequences) {
        out << ">1 1 " << weight << " 1 0\n";
        std::ranges::copy(sequence, ostream_iterator<Base>(out));
        out << "\n";
    }
}
