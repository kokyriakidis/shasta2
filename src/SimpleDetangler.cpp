#include "SimpleDetangler.hpp"
#include "Tangle.hpp"
#include "Tangle2.hpp"
#include "TangleMatrix2.hpp"
using namespace shasta;



bool SimpleDetangler::operator()(Tangle&)
{
    return false;
}



bool SimpleDetangler::operator()(Tangle2& tangle)
{
    const TangleMatrix2& tangleMatrix = *(tangle.tangleMatrix);

    // Check common coverage on all entrances and exits.
    for(const auto& entrance: tangleMatrix.entrances) {
        if(entrance.commonCoverage < minCommonCoverage) {
            return false;
        }
    }
    for(const auto& exit: tangleMatrix.exits) {
        if(exit.commonCoverage < minCommonCoverage) {
            return false;
        }
    }


    if(tangleMatrix.entrances.size() < 2) {
        return false;
    }
    if(tangleMatrix.exits.size() < 2) {
        return false;
    }

    // Check that each entrance will get a connection with at least one exit.
    for(uint64_t iEntrance=0; iEntrance<tangleMatrix.entrances.size(); iEntrance++) {
        bool isGood = false;
        for(uint64_t iExit=0; iExit<tangleMatrix.exits.size(); iExit++) {
            if(tangleMatrix.tangleMatrix[iEntrance][iExit].size() >= minDetangleCoverage) {
                isGood = true;
                break;
            }
        }
        if(not isGood) {
            return false;
        }
    }

    // Check that each exit will get a connection with at least one entrance.
    for(uint64_t iExit=0; iExit<tangleMatrix.exits.size(); iExit++) {
        bool isGood = false;
        for(uint64_t iEntrance=0; iEntrance<tangleMatrix.entrances.size(); iEntrance++) {
            if(tangleMatrix.tangleMatrix[iEntrance][iExit].size() >= minDetangleCoverage) {
                isGood = true;
                break;
            }
        }
        if(not isGood) {
            return false;
        }
    }


    // Do the detangling.
    // Connect the pairs with minDetangleCoverage.
    for(uint64_t iEntrance=0; iEntrance<tangleMatrix.entrances.size(); iEntrance++) {
        for(uint64_t iExit=0; iExit<tangleMatrix.exits.size(); iExit++) {
            if(tangleMatrix.tangleMatrix[iEntrance][iExit].size() >= minDetangleCoverage) {
                tangle.connect(iEntrance, iExit);
            }
        }

    }
    tangle.detangle();

    return true;

}



bool SimpleDetangler::operator()(Tangle3&)
{
    return false;
}
