// Shasta2.
#include "poastaWrapper.hpp"
#include "Base.hpp"
#include "SHASTA2_ASSERT.hpp"
using namespace shasta2;

// Poasta.
#include "poasta.h"

// Standard library.
#include "algorithm.hpp"
#include <cstring>



void shasta2::poasta(

    // The input sequences to be aligned.
    // They are presented to poasta in this order.
    const vector< vector<Base> >& sequences,

    // The consensus sequence and its coverage.
    vector< pair<Base, uint64_t> >& consensus,

    // The computed alignment.
    // Each element of the vector correspond to one of the input sequences,
    // in the same order.
    // These all have the same length, which equals the length of the aligned consensus.
    vector< vector<AlignedBase> >& alignment,

    // The aligned consensus.
    vector<AlignedBase>& alignedConsensus
)
{
    // Alignment parameters.
    const int mismatchScore = 4;
    const int gapOpenScore = 6;
    const int gapExtendScore = 2;

    // Create the PoastaGraph.
    PoastaGraph* graph = poasta_create_graph();

    // Add the sequences.
    string sequenceString;
    for(const vector<Base>& sequence: sequences) {
        sequenceString.clear();
        for(const Base base: sequence) {
            sequenceString.push_back(base.character());
        }
        poasta_add_sequence(graph, sequenceString.data(), sequence.size(),
            mismatchScore, gapOpenScore, gapExtendScore);
    }

    // Get the multiple sequence alignment.
    PoastaMsa msa = poasta_get_msa(graph);
    SHASTA2_ASSERT(msa.num_sequences == sequences.size());

    // Fill in the alignment.
    const uint64_t n = sequences.size();
    alignment.clear();
    for(uint64_t i=0; i<n; i++) {
        const char* sequenceCharacters = msa.sequences[i];
        vector<AlignedBase>& alignmentRow =alignment.emplace_back();
        for(uint64_t j=0; ; j++) {
            const char c = sequenceCharacters[j];
            if(c == 0) {
                break;
            }
            alignmentRow.push_back(AlignedBase::fromCharacter(c));
        }
    }



    // Get the alignment length and sanity check.
    const uint64_t alignmentLength = alignment.front().size();
    for(uint64_t i=1; i<alignment.size(); i++) {
        SHASTA2_ASSERT(alignment[i].size() == alignmentLength);
    }



    // Compute consensus and alignedConsensus.
    consensus.clear();
    alignedConsensus.resize(alignmentLength);
    // Loop over alignment positions.
    for(uint64_t i=0; i<alignmentLength; i++) {

        // Count bases at this position.
        array<uint64_t, 5> baseCount;
        fill(baseCount.begin(), baseCount.end(), 0);
        for(uint64_t j=0; j<n; j++) {
            const AlignedBase alignedBase = alignment[j][i];
            ++baseCount[alignedBase.value];
        }

        // Get the consensus AlignedBase at this position and its coverage.
        const auto it = std::max_element(baseCount.begin(), baseCount.end());
        const AlignedBase consensusBase = AlignedBase::fromInteger(uint64_t(it - baseCount.begin()));
        const uint64_t coverage = *it;

        // Store in alignedConsensus.
        alignedConsensus[i] = consensusBase;

        // Store in consensus.
        if(not consensusBase.isGap()) {
            consensus.push_back(make_pair(Base(consensusBase), coverage));
        }


    }

    // Free the poasta objects.
    poasta_free_msa(msa);
    poasta_free_graph(graph);
}



void shasta2::testPoasta1()
{
    // 1. Create a new graph
    PoastaGraph* graph = poasta_create_graph();

    // 2. Add sequences
    // Parameters: graph, sequence, length, mismatch_score, gap_open, gap_extend
    // Default scoring: mismatch=4, gap_open=6, gap_extend=2

    const char* seq1 = "ACGTACGT";
    poasta_add_sequence(graph, seq1, std::strlen(seq1), 4, 6, 2);

    const char* seq2 = "ACGTCGT";
    poasta_add_sequence(graph, seq2, std::strlen(seq2), 4, 6, 2);

    // 3. Get Multiple Sequence Alignment (MSA)
    PoastaMsa msa = poasta_get_msa(graph);

    std::cout << "MSA (" << msa.num_sequences << " sequences):" << std::endl;
    for (size_t i = 0; i < msa.num_sequences; ++i) {
        std::cout << msa.sequences[i] << std::endl;
    }

    // Free MSA memory
    poasta_free_msa(msa);

#if 0
    // 4. Get GFA output (optional)
    char* gfa = poasta_get_gfa(graph);
    if (gfa) {
        std::cout << "\nGFA Output:\n" << gfa << std::endl;
        free(gfa); // Use standard free() for the string
    }
#endif

    // 5. Free the graph
    poasta_free_graph(graph);


}


void shasta2::testPoasta2()
{
    // Define the sequences for the MSA.
    const uint64_t n = 3;
    const array<string, n> sequenceStrings = {"CTAGTT", "CTAAAGTGT", "ATAAAGTT"};
    vector< vector<Base> > sequences;
    for(const string& sequenceString: sequenceStrings) {
        vector<Base>& sequence = sequences.emplace_back();
        for(const char c: sequenceString) {
            sequence.emplace_back(Base::fromCharacter(c));
        }
    }

    // Write out the input sequences.
    cout << "Input sequences:" << endl;
    for(const vector<Base>& sequence: sequences) {
        for(const Base base: sequence) {
            cout << base;
        }
        cout << endl;
    }

    // Run the poasta wrapper.
    vector< pair<Base, uint64_t> > consensus;
    vector< vector<AlignedBase> > alignment;
    vector<AlignedBase> alignedConsensus;
    poasta(sequences, consensus, alignment, alignedConsensus);

    // Write out the alignment.
    cout << "Alignment:" << endl;
    for(const vector<AlignedBase>& alignmentRow: alignment) {
        for(const AlignedBase base: alignmentRow) {
            cout << base;
        }
        cout << endl;
    }

    // Write out the aligned consensus.
    cout << "Aligned consensus:" << endl;
    for(const AlignedBase& base: alignedConsensus) {
        cout << base;
    }
    cout << endl;

    // Write out the consensus.
    cout << "Consensus:" << endl;
    for(const auto& [base, coverage]: consensus) {
        cout << base;
    }
    cout << endl;
    for(const auto& [base, coverage]: consensus) {
        if(coverage < 10) {
            cout << coverage;
        } else {
            cout << "*";
        }
    }
    cout << endl;
}
