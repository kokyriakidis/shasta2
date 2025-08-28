// Shasta.
#include "LikelihoodRatioDetangler.hpp"
#include "GTest.hpp"
#include "Tangle1.hpp"
#include "TangleMatrix1.hpp"
using namespace shasta;



LikelihoodRatioDetangler::LikelihoodRatioDetangler(
    const double epsilon,
    const double maxLogP,
    const double minLogPDelta):
    epsilon(epsilon),
    maxLogP(maxLogP),
    minLogPDelta(minLogPDelta)
{}



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

        for(uint64_t i=0; i<entranceCount; i++) {
            cout << "Entrance " << tangle.assemblyGraph[tangle.entrances[i]].id << " oriented reads:" << endl;
            for(const TangleMatrix1::OrientedReadInfo& info: tangleMatrix.entranceOrientedReadInfos[i]) {
                cout << info.orientedReadId << " ";
            }
            cout << endl;
        }

        for(uint64_t j=0; j<exitCount; j++) {
            cout << "Exit " << tangle.assemblyGraph[tangle.exits[j]].id << " oriented reads:" << endl;
            for(const TangleMatrix1::OrientedReadInfo& info: tangleMatrix.exitOrientedReadInfos[j]) {
                cout << info.orientedReadId << " ";
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

