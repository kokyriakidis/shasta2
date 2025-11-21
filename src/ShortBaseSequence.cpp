#include "ShortBaseSequence.hpp"
#include "ShortBaseSequenceEdit.hpp"
#include "algorithm.hpp"
#include "SHASTA2_ASSERT.hpp"
using namespace shasta2;

#include <iomanip>



void shasta2::testShortBaseSequence()
{
    while(true) {

        // Read the k-mer string.
        string kmerString;
        cout << "Enter a k-mer." << endl;
        cin >> kmerString;
        const uint64_t k = kmerString.size();
        if(k > 128) {
            cout << "Can be at most 128 bases long." << endl;
        }

        // Create the k-mer.
        ShortBaseSequence128 kmer;
        for(uint64_t i=0; i<k; i++) {
            kmer.set(i, Base::fromCharacter(kmerString[i]));
        }

        kmer.write(cout, 128);
        cout << endl;

        for(uint64_t i=0; i<128; i++) {
            cout << kmer[i];
        }
        cout << endl;

#if 0
        // Test maxHomopolymerLength.
        cout << "maxHomopolymerLength returned " <<
            kmer.maxHomopolymerLength(k) << endl;

        // Test countExactRepeatCopies.
        cout << "countExactRepeatCopies<1> returned " <<
            kmer.countExactRepeatCopies<1>(k) << endl;
        cout << "countExactRepeatCopies<2> returned " <<
            kmer.countExactRepeatCopies<2>(k) << endl;
        cout << "countExactRepeatCopies<3> returned " <<
            kmer.countExactRepeatCopies<3>(k) << endl;
        cout << "countExactRepeatCopies<4> returned " <<
            kmer.countExactRepeatCopies<4>(k) << endl;
        cout << "countExactRepeatCopies<5> returned " <<
            kmer.countExactRepeatCopies<5>(k) << endl;
        cout << "countExactRepeatCopies<6> returned " <<
            kmer.countExactRepeatCopies<6>(k) << endl;
#endif
    }
}
