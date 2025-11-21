// Shasta2.
#include "poastaWrapper.hpp"
#include "Base.hpp"
#include "SHASTA2_ASSERT.hpp"
#include "simpleFastaRead.hpp"
#include "tmpDirectory.hpp"
using namespace shasta;

// Boost libraries.
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>

// Standard library.
#include <algorithm.hpp>
#include <filesystem>
#include "fstream.hpp"
#include <iterator.hpp>



void shasta::poasta(

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
    const uint64_t n = sequences.size();

    const string baseFileName = tmpDirectory() + to_string(boost::uuids::random_generator()());
    const string sequencesFileName = baseFileName + "-sequences.fasta";
    const string msaFileName = baseFileName + "-msa.fasta";

    // Write out the sequences to a fasta file.
    {
        ofstream fasta(sequencesFileName);
        for(uint64_t i=0; i<n; i++) {
            const vector<Base>& sequence = sequences[i];
            fasta << ">" << i << "\n";
            copy(sequence.begin(), sequence.end(), ostream_iterator<Base>(fasta));
            fasta << "\n";
        }
    }

    // Run poasta. There are no C/C++ bindings, so we have to invoke the executable.
    // This is only practical for long alignments.
    const string command = "poasta align -O fasta -o " + msaFileName + " " + sequencesFileName + " 2>/dev/null";
    const int commandStatus = std::system(command.c_str());
    SHASTA2_ASSERT(commandStatus == 0);
    std::filesystem::remove(sequencesFileName);

    // Read poasta output containing the alignment.
    {
        ifstream fasta(msaFileName);
        simpleFastaRead(fasta, alignment);
        SHASTA2_ASSERT(alignment.size() == sequences.size());
    }
    std::filesystem::remove(msaFileName);

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
}
