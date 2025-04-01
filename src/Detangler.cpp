#include "Detangler.hpp"
#include "Tangle.hpp"
using namespace shasta;



bool TrivialDetangler::operator()(
    Tangle& tangle)
{
    const TangleMatrix& tangleMatrix = tangle.tangleMatrix;

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


    if(tangleMatrix.entrances.size() != 2) {
        return false;
    }
    if(tangleMatrix.exits.size() != 2) {
        return false;
    }

    const bool isPositive00 = tangleMatrix.tangleMatrix[0][0] > 0;
    const bool isPositive01 = tangleMatrix.tangleMatrix[0][1] > 0;
    const bool isPositive10 = tangleMatrix.tangleMatrix[1][0] > 0;
    const bool isPositive11 = tangleMatrix.tangleMatrix[1][1] > 0;

    const bool isInPhase =
        isPositive00 and isPositive11 and (not isPositive01) and (not isPositive10);
    const bool isOutOfPhase =
        isPositive01 and isPositive10 and (not isPositive00) and (not isPositive11);

    if(not (isInPhase or isOutOfPhase)) {
        return false;
    }

    // Do the detangling.
    if(isInPhase) {
        tangle.connect(0, 0);
        tangle.connect(1, 1);
    } else if(isOutOfPhase) {
        tangle.connect(0, 1);
        tangle.connect(1, 0);
    }
    tangle.detangle();

    return true;

}
