#include "theseusWrapper.hpp"
#include "Base.hpp"
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



void shasta2::theseus(

    // The input sequences fixed on both sides, with their coverage.
    const vector< pair<vector<Base>, uint64_t> >& fixedSequences,

    // The input sequences fixed on the left only, with their coverage.
    const vector< pair<vector<Base>, uint64_t> >& leftFixedSequences,

    // The input sequences fixed on the right only, with their coverage.
    const vector< pair<vector<Base>, uint64_t> >& rightFixedSequences,

    // The consensus sequence and its coverage.
    vector<Base>& consensus,
    vector<AlignedBase>& alignedConsensus,
    vector<uint64_t>& coverage,

    // The alignment is only computed if computeAlignment is true.
    vector< vector<AlignedBase> >& alignment,
    bool computeAlignment
)
{


    // Use pericles default penalties.
    const int match = 0;
    const int mismatch = 2;
    const int gapo = 3;
    const int gape = 1;
    theseus::Penalties penalties(match, mismatch, gapo, gape);

    // Theseus heuristics.
    theseus::Heuristics heuristics;

    // Create the theseus aligner, passing in the first sequence fixed on both sides.
    SHASTA2_ASSERT(not fixedSequences.empty());
    const auto& [firstSequence, firstSequenceWeight] = fixedSequences.front();
    theseus::TheseusMSA aligner(penalties, heuristics,
        toString(firstSequence), int(firstSequenceWeight), false);

    // Pass to the aligner the remaining sequences fixed on both sides.
    for(uint64_t i=1; i<fixedSequences.size(); i++) {
        const auto& [sequence, weight] = fixedSequences[i];
        const bool densityDrop = false;
        const bool lagPruning = false;
        aligner.align(toString(sequence), int(weight), false, false, densityDrop, lagPruning);
    }

    // Pass to the aligner the sequences fixed on the left only
    for(const auto& [sequence, weight]:  leftFixedSequences) {
        const bool densityDrop = true;
        const bool lagPruning = false;
        aligner.align(toString(sequence), int(weight), false, true, densityDrop, lagPruning);
    }

    // Pass to the aligner the sequences fixed on the right only
    for(const auto& [sequence, weight]:  rightFixedSequences) {
        const bool densityDrop = true;
        const bool lagPruning = false;
        aligner.align(toString(sequence), int(weight), true, true, densityDrop, lagPruning);
    }

    // Get the consensus.
    vector<int> coverageInt;
    string consensusString;
    string alignedConsensusString;
    aligner.majority_voting_consensus(coverageInt, consensusString, alignedConsensusString);
    consensus = vectorOfBasesFromString(consensusString);
    alignedConsensus = vectorOfAlignedBasesFromString(alignedConsensusString);
    coverage.clear();
    std::ranges::copy(coverageInt, back_inserter(coverage));

    // Also compute the alignment, if requested.
    if(computeAlignment) {
        alignment.clear();
        std::ostringstream s;
        aligner.print_as_msa(s);
        const string& alignmentString = s.str();

        boost::tokenizer< boost::char_separator<char> > tokenizer(alignmentString, boost::char_separator<char>("\n"));

        // The alignment string has a header line plus an alignment line
        // for each of the input sequences.
        for(const string& line: tokenizer) {
            SHASTA2_ASSERT(not line.empty());
            if(line[0] != '>') {
                alignment.push_back(vectorOfAlignedBasesFromString(line));
                // cout << std::string_view(line) << "\n";
            }
        }

        // The last line contains aligned consensus.
        // We already got aligned consensus from the call
        // to majority_voting_consensus, so we cna ignore it.
        alignment.pop_back();
        SHASTA2_ASSERT(alignment.size() ==
            fixedSequences.size() + leftFixedSequences.size() + rightFixedSequences.size());
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




    vector<Base> consensus;
    vector<AlignedBase> alignedConsensus;
    vector<uint64_t> coverage;
    vector< vector<AlignedBase> > alignment;
    theseus(fixedSequences, leftFixedSequences, rightFixedSequences,
        consensus, alignedConsensus, coverage, alignment, true);

    cout << "Consensus with coverage:" << endl;
    std::ranges::copy(consensus, ostream_iterator<Base>(cout));
    cout << endl;
    for(const uint64_t c: coverage) {
        if(c < 10) {
            cout << c;
        } else if(c < 36) {
            cout << char(c -10 + 'A');
        } else {
            cout << "*";
        }
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
