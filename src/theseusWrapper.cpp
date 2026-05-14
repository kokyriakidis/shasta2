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

    // The consensus sequence.
    vector<Base>& consensus,

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

    // Use pericles default heuristics.
    const bool densityDrop = false;
    const bool lagPruning = false;
    theseus::Heuristics heuristics(densityDrop, lagPruning);

    // Create the theseus aligner, passing in the first sequence fixed on both sides.
    SHASTA2_ASSERT(not fixedSequences.empty());
    const auto& [firstSequence, firstSequenceWeight] = fixedSequences.front();
    theseus::TheseusMSA aligner(penalties, heuristics,
        toString(firstSequence), int(firstSequenceWeight), false);

    // Pass to the aligner the remaining sequences fixed on both sides.
    for(uint64_t i=1; i<fixedSequences.size(); i++) {
        const auto& [sequence, weight] = fixedSequences[i];
        aligner.align(toString(sequence), int(weight), false, false);
    }

    // Pass to the aligner the sequences fixed on the left only
    for(const auto& [sequence, weight]:  leftFixedSequences) {
        aligner.align(toString(sequence), int(weight), false, true);
    }

    // Pass to the aligner the sequences fixed on the right only
    for(const auto& [sequence, weight]:  rightFixedSequences) {
        aligner.align(toString(sequence), int(weight), true, true);
    }

    // Get the consensus.
    consensus = vectorOfBasesFromString(aligner.get_majority_voting_consensus_sequence());


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
        SHASTA2_ASSERT(alignment.size() ==
            fixedSequences.size() + leftFixedSequences.size() + rightFixedSequences.size());
    }
}



void shasta2::testTheseus()
{
    const vector< pair<vector<Base>, uint64_t> > fixedSequences =
    {
        {vectorOfBasesFromString("TAGGGATTGATAAAAGGCTTTCCAGAAGA"), 5},
        {vectorOfBasesFromString("TAGGGATTGATAAAGGCTTTCCAGAAGA"), 4},
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
    vector< vector<AlignedBase> > alignment;
    theseus(fixedSequences, leftFixedSequences, rightFixedSequences, consensus, alignment, true);

    cout << "Consensus is:" << endl;
    std::ranges::copy(consensus, ostream_iterator<Base>(cout));
    cout << endl;


    cout << "Alignment:" << endl;
    for(const vector<AlignedBase>& alignmentRow: alignment) {
        for(const AlignedBase base: alignmentRow) {
            cout << base;
        }
        cout << endl;
    }
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
        out << ">0 0 " << weight << "\n";
        std::ranges::copy(sequence, ostream_iterator<Base>(out));
        out << "\n";
    }

    for(const auto& [sequence, weight]: leftFixedSequences) {
        out << ">0 1 " << weight << "\n";
        std::ranges::copy(sequence, ostream_iterator<Base>(out));
        out << "\n";
    }

    for(const auto& [sequence, weight]: rightFixedSequences) {
        out << ">1 1 " << weight << "\n";
        std::ranges::copy(sequence, ostream_iterator<Base>(out));
        out << "\n";
    }
}
