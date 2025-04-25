#include "simpleFastaRead.hpp"
#include "Base.hpp"
#include "SHASTA_ASSERT.hpp"
using namespace shasta;

#include <iostream.hpp>
#include <string.hpp>



void shasta::simpleFastaRead(istream& fasta, vector< vector<AlignedBase> >& sequences)
{
    sequences.clear();

    string line;

    // Read one line at a time.
    while(fasta.good()) {
        getline(fasta, line);

        // If the line is empty, ignore it.
        if(line.empty()) {
            continue;
        }

        // If the line begins with ">", start a new sequence
        // and ignore the rest of the line.
        if(line[0] == '>') {
            sequences.resize(sequences.size() + 1);
            continue;
        }

        // In all other cases, append the rest of the line to the
        // last sequence.
        SHASTA_ASSERT(not sequences.empty());
        vector<AlignedBase>& sequence = sequences.back();
        for(const char c: line) {
            sequence.push_back(AlignedBase::fromCharacter(c));
        }
    }
}
