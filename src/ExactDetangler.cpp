#include "ExactDetangler.hpp"
#include "Tangle3.hpp"
#include "TangleMatrix3.hpp"
using namespace shasta;



bool ExactDetangler::operator()(Tangle3& tangle)
{
    const TangleMatrix3& tangleMatrix = *(tangle.tangleMatrix);

    for(uint64_t iEntrance=0; iEntrance<tangleMatrix.entrances.size(); iEntrance++) {
        for(uint64_t iExit=0; iExit<tangleMatrix.exits.size(); iExit++) {
            const uint64_t coverage = tangleMatrix.tangleMatrix[iEntrance][iExit].size();
            if(coverage > 0) {
            	tangle.connect(iEntrance, iExit);
            }
        }
    }

    tangle.detangle();
    return true;
}
