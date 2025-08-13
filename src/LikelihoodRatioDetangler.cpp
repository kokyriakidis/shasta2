// Shasta.
#include "LikelihoodRatioDetangler.hpp"
#include "GTest.hpp"
#include "Tangle.hpp"
#include "Tangle1.hpp"
#include "TangleMatrix.hpp"
#include "TangleMatrix1.hpp"
using namespace shasta;



LikelihoodRatioDetangler::LikelihoodRatioDetangler(
    uint64_t minCommonCoverage,
    const double epsilon,
    const double maxLogP,
    const double minLogPDelta,
    uint64_t detangleHighCoverageThreshold,
    bool useExtendedTangleMatrix):
    minCommonCoverage(minCommonCoverage),
    epsilon(epsilon),
    maxLogP(maxLogP),
    minLogPDelta(minLogPDelta),
    detangleHighCoverageThreshold(detangleHighCoverageThreshold),
    useExtendedTangleMatrix(useExtendedTangleMatrix)
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

    vector< vector<uint64_t> > tangleMatrixCoverage;
    tangleMatrix.getTangleMatrixCoverage(tangleMatrixCoverage);



    // Run the likelihood ratio test.
    shared_ptr<GTest> gTestPointer;
    if(useExtendedTangleMatrix) {
        vector< vector<double> > extendedTangleMatrix;
        tangle.computeExtendedTangleMatrix(extendedTangleMatrix);
        gTestPointer = make_shared<GTest>(extendedTangleMatrix, epsilon);
    } else {
        gTestPointer = make_shared<GTest>(tangleMatrixCoverage, epsilon);
    }
    const GTest& gTest = *gTestPointer;
    if(not gTest.success) {
        if(debug) {
            cout << "Not detangling because the G-test failed." << endl;
        }
        return false;
    }



    // If the best G is not good enough, do nothing.
    const double bestG = gTest.hypotheses.front().G;
    if(bestG > maxLogP) {
        if(debug) {
            cout << "Not detangling because the best G is too large." << endl;
        }
        return false;
    }

    // Also check the second best G.
    if(gTest.hypotheses.size() >= 2) {
        const double secondBestG = gTest.hypotheses[1].G;
        const double GDelta = secondBestG - bestG;
        if(GDelta < minLogPDelta) {
            if(debug) {
                cout << "Not detangling because the second best G is too small." << endl;
            }
            return false;
        }
    }
    const vector< vector<bool> >& bestConnectivityMatrix = gTest.hypotheses.front().connectivityMatrix;



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

    // Figure out if the best connectivity matrix is backward injective (only one entrance for each exit).
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

    // If detangling would generate a connection with coverage
    // less than uint64_t detangleHighCoverageThreshold, don't detangle.
    for(uint64_t iEntrance=0; iEntrance<entranceCount; iEntrance++) {
        for(uint64_t iExit=0; iExit<exitCount; iExit++) {
            if(bestConnectivityMatrix[iEntrance][iExit]) {
                if(tangleMatrixCoverage[iEntrance][iExit] < detangleHighCoverageThreshold) {
                    if(debug) {
                        cout << "Not detangling to avoid generating an assembly step "
                            "with coverage " << tangleMatrixCoverage[iEntrance][iExit] << endl;
                    }
                    return false;
                }
            }
        }
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



bool LikelihoodRatioDetangler::operator()(Tangle1& tangle)
{
    using edge_descriptor = AssemblyGraph::edge_descriptor;

    const uint64_t entranceCount = tangle.entrances.size();
    const uint64_t exitCount = tangle.exits.size();
    const TangleMatrix1& tangleMatrix = tangle.tangleMatrix();

    if(debug) {
        cout << "Detangling begins for a tangle with " << entranceCount <<
            " entrances and " << exitCount << " exits." << endl;

        cout << "Entrances:";
        for(const edge_descriptor e: tangle.entrances) {
            cout << " " << tangle.assemblyGraph[e].id;
        }
        cout << endl;

        cout << "Exits:";
        for(const edge_descriptor e: tangle.exits) {
            cout << " " << tangle.assemblyGraph[e].id;
        }
        cout << endl;

        cout << "Tangle matrix:" << endl;
        for(uint64_t i=0; i<entranceCount; i++) {
            for(uint64_t j=0; j<exitCount; j++) {
                cout << tangleMatrix.tangleMatrix[i][j] << " ";
            }
            cout << endl;
        }
    }

    // Run the likelihood ratio test.
    const GTest gTest(tangleMatrix.tangleMatrix, epsilon);
    if(not gTest.success) {
        if(debug) {
            cout << "Not detangling because the G-test failed." << endl;
        }
        return false;
    }



    // If the best G is not good enough, do nothing.
    const double bestG = gTest.hypotheses.front().G;
    if(debug) {
        cout << "Best G " << bestG << endl;
    }
    if(bestG > maxLogP) {
        if(debug) {
            cout << "Not detangling because the best G is too large." << endl;
        }
        return false;
    }

    // Also check the second best G.
    if(gTest.hypotheses.size() >= 2) {
        const double secondBestG = gTest.hypotheses[1].G;
        if(debug) {
            cout << "Second best G " << secondBestG << endl;
        }
        const double GDelta = secondBestG - bestG;
        if(GDelta < minLogPDelta) {
            if(debug) {
                cout << "Not detangling because the second best G is too small." << endl;
            }
            return false;
        }
    }
    const vector< vector<bool> >& bestConnectivityMatrix = gTest.hypotheses.front().connectivityMatrix;



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

    // Figure out if the best connectivity matrix is backward injective (only one entrance for each exit).
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



    // If getting here, we will detangle using the connectivity matrix for the best hypothesis.
    if(debug) {
        cout << "Connectivity matrix of best hypothesis:" << endl;
        for(uint64_t i=0; i<entranceCount; i++) {
            for(uint64_t j=0; j<exitCount; j++) {
                cout << uint64_t(bestConnectivityMatrix[i][j]) << " ";
            }
            cout << endl;
        }

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

