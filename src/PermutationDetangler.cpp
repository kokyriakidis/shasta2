#include "PermutationDetangler.hpp"
#include "Tangle.hpp"
using namespace shasta;



bool PermutationDetangler::operator()(
    Tangle& tangle)
{
    const TangleMatrix& tangleMatrix = tangle.tangleMatrix;
    const vector< vector<uint64_t> >& matrix = tangleMatrix.tangleMatrix;

    // The PermutationDetangler only works on Tangles where the
    // numbers of Entrances and Exits are the same.
    const uint64_t n = tangleMatrix.entrances.size();
    if(tangleMatrix.exits.size() != n) {
        return false;
    }

    // Check common coverage on all entrances and exits.
    for(const TangleMatrix::Entrance& entrance: tangleMatrix.entrances) {
        if(entrance.commonCoverage < minCommonCoverage) {
            return false;
        }
    }
    for(const TangleMatrix::Exit& exit: tangleMatrix.exits) {
        if(exit.commonCoverage < minCommonCoverage) {
            return false;
        }
    }

    if(debug) {
        cout << "PermutationDetangler called with the following TangleMatrix:" << endl;
        for(uint64_t iEntrance=0; iEntrance<n; iEntrance++) {
            for(uint64_t iExit=0; iExit<n; iExit++) {
                cout << iEntrance << " " << iExit << " " << matrix[iEntrance][iExit] << endl;
            }
        }
    }



    return false;

}
