// Shasta.
#include "LikelihoodRatioDetangler.hpp"
// #include "Base.hpp"
// #include "orderPairs.hpp"
#include "Tangle.hpp"
#include "TangleMatrix.hpp"
using namespace shasta;



LikelihoodRatioDetangler::LikelihoodRatioDetangler(
    uint64_t minCommonCoverage,
    const double epsilon,
    const double maxLogP,
    const double minLogPDelta):
    minCommonCoverage(minCommonCoverage),
    epsilon(epsilon),
    maxLogP(maxLogP),
    minLogPDelta(minLogPDelta)
{}



bool LikelihoodRatioDetangler::operator()(Tangle& tangle, bool doDetangle)
{
    TangleMatrix& tangleMatrix = *(tangle.tangleMatrix);
    const uint64_t entranceCount = tangleMatrix.entrances.size();
    const uint64_t exitCount = tangleMatrix.exits.size();

    if((entranceCount < 2) or (exitCount < 2)) {
        if(debug) {
            cout << "Not detangling due to in-degree or out-degree." << endl;
        }
        return false;
    }

    if(debug) {
        cout << "Working on a tangle with " << entranceCount << " entrances and " << exitCount << " exits." << endl;
        cout << "Entrances:";
        for(const auto& entrance: tangleMatrix.entrances) {
            cout << " " << tangle.assemblyGraph[entrance.e].id;
        }
        cout << endl;
        cout << "Exits:";
        for(const auto& exit: tangleMatrix.exits) {
            cout << " " << tangle.assemblyGraph[exit.e].id;
        }
        cout << endl;
        cout << "Tangle matrix:" << endl;
        for(uint64_t iEntrance=0; iEntrance<entranceCount; iEntrance++) {
            for(uint64_t iExit=0; iExit<exitCount; iExit++) {
                const uint64_t coverage = tangleMatrix.tangleMatrix[iEntrance][iExit].size();
                cout << coverage << " ";
            }
            cout << endl;
        }
    }

    // Check common coverage on all entrances and exits.
    for(const auto& entrance: tangleMatrix.entrances) {
        if(entrance.commonCoverage < minCommonCoverage) {
            if(debug) {
                cout << "Insufficient common coverage." << endl;
            }
            return false;
        }
    }
    for(const auto& exit: tangleMatrix.exits) {
        if(exit.commonCoverage < minCommonCoverage) {
            if(debug) {
                cout << "Insufficient common coverage." << endl;
            }
            return false;
        }
    }

    // Run the likelihood ratio test.
    if(not tangleMatrix.gTest(epsilon)) {
        if(debug) {
            cout << "Not detangling because the tangle matrix is too big." << endl;
        }
        return false;
    }




    // If the best G is not good enough, do nothing.
    const double bestG = tangleMatrix.hypotheses.front().G;
    if(bestG > maxLogP) {
        if(debug) {
            cout << "Not detangling because the best G is too large." << endl;
        }
        return false;
    }

    // Also check the second best G.
    if(tangleMatrix.hypotheses.size() >= 2) {
        const double secondBestG = tangleMatrix.hypotheses[1].G;
        const double GDelta = secondBestG - bestG;
        if(GDelta < minLogPDelta) {
            if(debug) {
                cout << "Not detangling because the second best G is too small." << endl;
            }
            return false;
        }
    }

    // Store the connect pairs and, if requested, detangle.
    vector< vector<bool> >& bestConnectivityMatrix = tangleMatrix.hypotheses.front().connectivityMatrix;
    for(uint64_t iEntrance=0; iEntrance<entranceCount; iEntrance++) {
        for(uint64_t iExit=0; iExit<exitCount; iExit++) {
            if(bestConnectivityMatrix[iEntrance][iExit]) {
                tangle.connect(iEntrance, iExit);
            }
        }
    }
    if(doDetangle) {
        tangle.detangle();
    }

    return true;
}
