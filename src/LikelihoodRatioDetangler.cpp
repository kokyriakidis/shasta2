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



bool LikelihoodRatioDetangler::operator()(Tangle& tangle)
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

#if 0
        if( (tangleMatrix.entrances.size() == 2) and
            (tangleMatrix.exits.size() == 2) and
            (tangle.assemblyGraph[tangleMatrix.entrances[0].e].id == 33053) and
            (tangle.assemblyGraph[tangleMatrix.entrances[1].e].id == 98033) and
            (tangle.assemblyGraph[tangleMatrix.exits    [0].e].id == 13086) and
            (tangle.assemblyGraph[tangleMatrix.exits    [1].e].id == 97642)) {
            tangle.assemblyGraph.write("Y");
        }
#endif

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
    const vector< vector<bool> >& bestConnectivityMatrix = tangleMatrix.hypotheses.front().connectivityMatrix;



    // Figure out if the best connectivity matrix is forward injective (only one exit for each entrance).
    bool isForwardInjective = true;
    for(uint64_t iEntrance=0; iEntrance<entranceCount; iEntrance++) {
        uint64_t count = 0;
        for(uint64_t iExit=0; iExit<exitCount; iExit++) {
            if(bestConnectivityMatrix[iEntrance][iExit]) {
                ++count;
            }
        }
        if(count > 1) {
            isForwardInjective = false;
            break;
        }
    }

    // Figure out if the best connectivity matrix is backward injective (only one enytrance for each exit).
    bool isBackwardInjective = true;
    for(uint64_t iExit=0; iExit<exitCount; iExit++) {
        uint64_t count = 0;
        for(uint64_t iEntrance=0; iEntrance<entranceCount; iEntrance++) {
            if(bestConnectivityMatrix[iEntrance][iExit]) {
                ++count;
            }
        }
        if(count > 1) {
            isBackwardInjective = false;
            break;
        }
    }

    if(not(isForwardInjective or isBackwardInjective)) {
        if(debug) {
            cout << "Not detangling because the best connectivity matrix is not forward or backward injective." << endl;
        }
        return false;
    }



    // Store the connect pairs and detangle.
    for(uint64_t iEntrance=0; iEntrance<entranceCount; iEntrance++) {
        for(uint64_t iExit=0; iExit<exitCount; iExit++) {
            if(bestConnectivityMatrix[iEntrance][iExit]) {
                tangle.connect(iEntrance, iExit);
                if(debug) {
                    cout << "Connecting for detangling: " << iEntrance << " " << iExit << endl;
                }
            }
        }
    }
    tangle.detangle();

    return true;
}
