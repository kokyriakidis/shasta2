#include "abpoaWrapper.hpp"
#include "Base.hpp"
#include "SHASTA2_ASSERT.hpp"
using namespace shasta2;

#include "abpoa.h"


// C++ wrapper to abpoa.
void shasta2::abpoa(

    // The input sequences to be aligned.
    // They are presented to abpoa in this order.
    const vector< vector<Base> >& sequences,

    // The consensus sequence and its coverage.
    vector< pair<Base, uint64_t> >& consensus,

    // The computed alignment.
    // Each element of the vector correspond to one of the input sequences,
    // in the same order.
    // These all have the same length, which equals the length of the aligned consensus.
    vector< vector<AlignedBase> >& alignment,

    // The aligned consensus.
    vector<AlignedBase>& alignedConsensus,

    // Consensus is always computed.
    // Alignment and alignedConsensus are only computed if this set to true.
    bool computeAlignment
)
{
    // Set up abpoa.
    abpoa_t* ab = abpoa_init();
    abpoa_para_t* abpt = abpoa_init_para();
    // alignment parameters
    // abpt->align_mode = 0; // 0:global 1:local, 2:extension
    // abpt->mat_fn = strdup("HOXD70.mtx"); abpt->use_score_matrix = 1; // score matrix instead of constant match/mismatch score
    // abpt->match = 2;      // match score
    // abpt->mismatch = 4;   // mismatch penalty
    // abpt->gap_mode = ABPOA_CONVEX_GAP; // gap penalty mode
    // abpt->gap_open1 = 4;  // gap open penalty #1
    // abpt->gap_ext1 = 2;   // gap extension penalty #1
    // abpt->gap_open2 = 24; // gap open penalty #2
    // abpt->gap_ext2 = 1;   // gap extension penalty #2
                             // gap_penalty = min{gap_open1 + gap_len * gap_ext1, gap_open2 + gap_len * gap_ext2}
    // abpt->bw = 10;        // extra band used in adaptive banded DP
    // abpt->bf = 0.01;

    // output options
    if(computeAlignment) {
        abpt->out_msa = 1; // generate Row-Column multiple sequence alignment(RC-MSA), set 0 to disable
    } else {
        abpt->out_msa = 0; // generate Row-Column multiple sequence alignment(RC-MSA), set 0 to disable
    }
    abpt->out_cons = 1; // generate consensus sequence, set 0 to disable
    abpt->w = 6, abpt->k = 9; abpt->min_w = 10; // minimizer-based seeding and partition
    abpt->progressive_poa = 0;
    abpt->max_n_cons = 1; // to generate 1 consensus sequences
    // abpt->sub_aln = 1;

    // Use majority voting to compute consensus. The default is to use the MSA graph,
    // and that sometimes generates unexpected consensuses.
    // See here for some discussion https://github.com/yangao07/abPOA/issues/67.
    abpt->cons_algrm = 1;

    abpoa_post_set_para(abpt);
    abpt->use_qv = 1;



    // The number of input sequences.
    const uint64_t n = sequences.size();

    // Vector to contain pointers to the input sequences.
    // This will require a reinterpret_cast from Base to uint8_t.
    vector<uint8_t*> sequencePointers;

    // Vector to contain the lengths of the input sequences.
    vector<int> sequenceLengths;

    for(const vector<Base>& sequence: sequences) {
        sequencePointers.push_back(reinterpret_cast<uint8_t*>(const_cast<Base*>(&sequence[0])));
        sequenceLengths.push_back(int(sequence.size()));
    }
    int* sequenceLengthsPointer = &sequenceLengths[0];



    // Call abpoa.
    abpoa_msa(ab, abpt, int(n), 0, sequenceLengthsPointer, &sequencePointers[0], 0, 0);

    // Store the consensus.
    abpoa_cons_t* abc = ab->abc;
    const int consensusLength = *abc->cons_len;
    uint8_t* consensusPointer = *(abc->cons_base);
    int* consensusCoverage = *(abc->cons_cov);
    consensus.clear();
    for(int i=0; i<consensusLength; i++) {
        consensus.push_back(make_pair(Base::fromInteger(consensusPointer[i]), consensusCoverage[i]));
    }

    // Store the alignment.
    const uint64_t alignmentLength = abc->msa_len;
    alignment.clear();
    alignment.resize(n, vector<AlignedBase>(alignmentLength));
    for(uint64_t i=0; i<n; i++) {
        for(uint64_t j=0; j<alignmentLength; j++) {
            const uint8_t b = abc->msa_base[i][j];
            AlignedBase base = (b < 4) ? AlignedBase::fromInteger(b) : AlignedBase::gap();
            alignment[i][j] = base;
        }

    }


    // Store aligned consensus.
    alignedConsensus.resize(alignmentLength);
    for(uint64_t j=0; j<alignmentLength; j++) {
        const uint8_t b = abc->msa_base[n][j];
        AlignedBase base = (b < 4) ? AlignedBase::fromInteger(b) : AlignedBase::gap();
        alignedConsensus[j] = base;
    }

    abpoa_free(ab);
    abpoa_free_para(abpt);
}
