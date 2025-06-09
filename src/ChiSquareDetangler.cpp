// Shasta.
#include "ChiSquareDetangler.hpp"
#include "Base.hpp"
#include "orderPairs.hpp"
#include "Tangle.hpp"
#include "TangleMatrix.hpp"
using namespace shasta;

// Boost libraries.
#include <boost/math/distributions/chi_squared.hpp>



ChiSquareDetangler::ChiSquareDetangler(
    uint64_t minCommonCoverage,
    const double epsilon,
    const double maxLogP,
    const double minLogPDelta):
    minCommonCoverage(minCommonCoverage),
    epsilon(epsilon),
    maxLogP(maxLogP),
    minLogPDelta(minLogPDelta)
{}



bool ChiSquareDetangler::operator()(Tangle& tangle)
{
    const TangleMatrix& tangleMatrix = *(tangle.tangleMatrix);
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



    // Special case : 2 by 2 tangle matrix.
    if((entranceCount == 2) and (exitCount == 2)) {

        // Special case: 2 by 2 tangle matrix, diagonal.
        if((tangleMatrix.tangleMatrix[0][1].size() == 0) and (tangleMatrix.tangleMatrix[1][0].size() == 0)) {
            tangle.connect(0, 0);
            tangle.connect(1, 1);
            tangle.detangle();
            if(debug) {
                cout << "Special case: 2 by 2, in phase." << endl;
            }
            return true;
        }

        // Special case: 2 by 2 tangle matrix, off-diagonal.
        if((tangleMatrix.tangleMatrix[0][0].size() == 0) and (tangleMatrix.tangleMatrix[1][1].size() == 0)) {
            tangle.connect(0, 1);
            tangle.connect(1, 0);
            tangle.detangle();
            if(debug) {
                cout << "Special case: 2 by 2, out of phase." << endl;
            }
            return true;
        }
    }



    // Code for the general case follows.



    // For now limit this to 2 by 2, 2 by 3, and 3 by 3.
    const uint64_t totalTangleMatrixEntryCount = entranceCount * exitCount;
    if(totalTangleMatrixEntryCount > 9) {
        if(debug) {
            cout << "Not detangling because the tangle matrix is too big." << endl;
        }
        return false;
    }



    // Chi squared distribution used below to test each permutation.
    boost::math::chi_squared_distribution chi2Distribution(double(totalTangleMatrixEntryCount - 1));



    // Compute total common coverage.
    uint64_t totalCommonCoverageOnEntrances = 0;
    for(const auto& entrance: tangleMatrix.entrances) {
        totalCommonCoverageOnEntrances += entrance.commonCoverage;
    }
    uint64_t totalCommonCoverageOnExits = 0;
    for(const auto& exit: tangleMatrix.exits) {
        totalCommonCoverageOnExits += exit.commonCoverage;
    }
    SHASTA_ASSERT(totalCommonCoverageOnEntrances == totalCommonCoverageOnExits);
    const uint64_t totalCommonCoverage = totalCommonCoverageOnEntrances;



    // Compute what the tangle matrix would be under entirely random assumptions.
    vector< vector<double> > randomTangleMatrix(entranceCount, vector<double>(exitCount, 0.));
    for(uint64_t iEntrance=0; iEntrance<entranceCount; iEntrance++) {
        for(uint64_t iExit=0; iExit<exitCount; iExit++) {
            randomTangleMatrix[iEntrance][iExit] =
                double(tangleMatrix.entrances[iEntrance].commonCoverage) *
                double(tangleMatrix.exits[iExit].commonCoverage) /
                double(totalCommonCoverage);
        }
    }

    if(debug) {
        cout << "Tangle matrix under random assumptions:" << endl;
        for(uint64_t iEntrance=0; iEntrance<entranceCount; iEntrance++) {
            for(uint64_t iExit=0; iExit<exitCount; iExit++) {
                cout << randomTangleMatrix[iEntrance][iExit] << " ";
            }
            cout << endl;
        }
    }

    // Some matrices used in the loop below.
    vector< vector<double> > idealTangleMatrixA(entranceCount, vector<double>(exitCount));
    vector< vector<double> > idealTangleMatrixB(entranceCount, vector<double>(exitCount));
    vector< vector<double> > idealTangleMatrix(entranceCount, vector<double>(exitCount));
    vector< vector<double> > expectedTangleMatrix(entranceCount, vector<double>(exitCount));

    // Keep track of the logP for the conectivity matrices for which it was not infinity.
    vector<pair<double, vector< vector<bool> > > > table;


    // The connectivity matrix contains true for entrance/exit pairs to be connected.
    // Try all N possible connectivity matrices and keep the one that gives the best chi squared.
    const uint64_t N = 1ULL << totalTangleMatrixEntryCount;
    vector< vector<bool> > connectivityMatrix(entranceCount, vector<bool>(exitCount));
    if(debug) {
        cout << "General case. Trying " << N << " possible connectivity matrices." << endl;
    }
    for(uint64_t connectivityInteger=0; connectivityInteger<N; connectivityInteger++) {

        // Use the bits of connectivityInteger to construct the connectivity matrix.
        uint64_t mask = 1;
        for(uint64_t iEntrance=0; iEntrance<entranceCount; iEntrance++) {
            for(uint64_t iExit=0; iExit<exitCount; iExit++) {
                connectivityMatrix[iEntrance][iExit] = ((connectivityInteger & mask) != 0);
                mask = mask << 1;
            }
        }

        if(false) {
            cout << "Trying connectivity matrix:" << endl;
            for(uint64_t iEntrance=0; iEntrance<entranceCount; iEntrance++) {
                for(uint64_t iExit=0; iExit<exitCount; iExit++) {
                    cout << connectivityMatrix[iEntrance][iExit] << " ";
                }
                cout << endl;
            }
        }



        // Compute the tangle matrix we would see assuming this connectivity matrix
        // and no errors. This can be done in two ways:
        // - Equally distributing common coverage at each entrance among all the
        //   exits for which the connectivity matrix is true for that entrance (idealTangleMatrixA).
        // - Equally distributing common coverage at each exit among all the
        //   entrances for which the connectivity matrix is true for that entrance (idealTangleMatrixB).
        // We average the two ways to compute the idealTangleMatrix.

        // Compute idealTangleMatrixA.
        for(uint64_t iEntrance=0; iEntrance<entranceCount; iEntrance++) {
            uint64_t nonZeroCount = 0;
            for(uint64_t iExit=0; iExit<exitCount; iExit++) {
                if(connectivityMatrix[iEntrance][iExit]) {
                    ++nonZeroCount;
                }
            }
            if(nonZeroCount == 0) {
                for(uint64_t iExit=0; iExit<exitCount; iExit++) {
                    idealTangleMatrixA[iEntrance][iExit] = 0.;
                }
            } else {
                const double value = double(tangleMatrix.entrances[iEntrance].commonCoverage) / double(nonZeroCount);
                for(uint64_t iExit=0; iExit<exitCount; iExit++) {
                    if(connectivityMatrix[iEntrance][iExit]) {
                        idealTangleMatrixA[iEntrance][iExit] = value;
                    } else {
                        idealTangleMatrixA[iEntrance][iExit] = 0.;
                    }
                }
            }
        }

        if(false) {
            cout << "Ideal tangle matrix computed using entrances:" << endl;
            for(uint64_t iEntrance=0; iEntrance<entranceCount; iEntrance++) {
                for(uint64_t iExit=0; iExit<exitCount; iExit++) {
                    cout << idealTangleMatrixA[iEntrance][iExit] << " ";
                }
                cout << endl;
            }
        }

        // Compute idealTangleMatrixB.
        for(uint64_t iExit=0; iExit<exitCount; iExit++) {
            uint64_t nonZeroCount = 0;
            for(uint64_t iEntrance=0; iEntrance<entranceCount; iEntrance++) {
                if(connectivityMatrix[iEntrance][iExit]) {
                    ++nonZeroCount;
                }
            }
            if(nonZeroCount == 0) {
                for(uint64_t iEntrance=0; iEntrance<entranceCount; iEntrance++) {
                    idealTangleMatrixB[iEntrance][iExit] = 0.;
                }
            } else {
                const double value = double(tangleMatrix.exits[iExit].commonCoverage) / double(nonZeroCount);
                for(uint64_t iEntrance=0; iEntrance<entranceCount; iEntrance++) {
                    if(connectivityMatrix[iEntrance][iExit]) {
                        idealTangleMatrixB[iEntrance][iExit] = value;
                    } else {
                        idealTangleMatrixB[iEntrance][iExit] = 0.;
                    }
                }
            }
        }

        if(false) {
            cout << "Ideal tangle matrix computed using exits:" << endl;
            for(uint64_t iEntrance=0; iEntrance<entranceCount; iEntrance++) {
                for(uint64_t iExit=0; iExit<exitCount; iExit++) {
                    cout << idealTangleMatrixB[iEntrance][iExit] << " ";
                }
                cout << endl;
            }
        }


        // Compute the ideal tangle matrix by averaging idealTangleMatrixA and idealTangleMatrixB.
        for(uint64_t iEntrance=0; iEntrance<entranceCount; iEntrance++) {
            for(uint64_t iExit=0; iExit<exitCount; iExit++) {
                idealTangleMatrix[iEntrance][iExit] =  0.5 *
                    (idealTangleMatrixA[iEntrance][iExit] + idealTangleMatrixB[iEntrance][iExit]);
            }
        }

        if(false) {
            cout << "Ideal tangle matrix:" << endl;
            for(uint64_t iEntrance=0; iEntrance<entranceCount; iEntrance++) {
                for(uint64_t iExit=0; iExit<exitCount; iExit++) {
                    cout << idealTangleMatrix[iEntrance][iExit] << " ";
                }
                cout << endl;
            }
        }


        // Compute the expected tangle matrix given this connectivity matrix.
        for(uint64_t iEntrance=0; iEntrance<entranceCount; iEntrance++) {
            for(uint64_t iExit=0; iExit<exitCount; iExit++) {
                expectedTangleMatrix[iEntrance][iExit] =
                    (1. - epsilon) * idealTangleMatrix [iEntrance][iExit] +
                    epsilon        * randomTangleMatrix[iEntrance][iExit];
            }
        }

        if(false) {
            cout << "Expected non-ideal tangle matrix for this connectivity matrix:" << endl;
            for(uint64_t iEntrance=0; iEntrance<entranceCount; iEntrance++) {
                for(uint64_t iExit=0; iExit<exitCount; iExit++) {
                    cout << expectedTangleMatrix[iEntrance][iExit] << " ";
                }
                cout << endl;
            }
        }
        // Do a chi-square test f the observed tangle matrix against the expected one.
        double chi2 = 0.;
        for(uint64_t iEntrance=0; iEntrance<entranceCount; iEntrance++) {
            for(uint64_t iExit=0; iExit<exitCount; iExit++) {
                const double expectedCoverage = expectedTangleMatrix[iEntrance][iExit];
                const double actualCoverage = double(tangleMatrix.tangleMatrix[iEntrance][iExit].size());
                const double delta = actualCoverage - expectedCoverage;
                chi2 += delta * delta / expectedCoverage;
            }
        }
        const double cumulativeChi2 = cdf(complement(chi2Distribution, chi2));
        if(cumulativeChi2 == 0.) {
            if(false) {
                cout << "Chi square cumulative is zero." << endl;
            }
        }

        if(cumulativeChi2 > 0.) {
            // Compute the corresponding logP in decibels.
            const double logP = - 10. * log10(cdf(complement(chi2Distribution, chi2)));
            table.push_back(make_pair(logP, connectivityMatrix));

            if(false) {
                cout << "Chi square for this connectivity matrix: " << chi2 <<
                    ", logP " << logP << " dB." << endl;
            }
        }

    // End of loop over connectivity matrices.
    }


    // Sort our table by increasing logP.
    sort(table.begin(), table.end(), OrderPairsByFirstOnly<double, vector< vector<bool> > >());

    if(debug) {
        if(table.empty()) {
            cout << "No viable connectivity matrices were found." << endl;
        } else {

            cout << "Best connectivity matrix has logP = " << table.front().first << " dB:" << endl;
            for(uint64_t iEntrance=0; iEntrance<entranceCount; iEntrance++) {
                for(uint64_t iExit=0; iExit<exitCount; iExit++) {
                    cout << table.front().second[iEntrance][iExit] << " ";
                }
                cout << endl;
            }

            if(table.size() > 1) {
                cout << "Second best connectivity matrix has logP = " << table[1].first << " dB:" << endl;
                for(uint64_t iEntrance=0; iEntrance<entranceCount; iEntrance++) {
                    for(uint64_t iExit=0; iExit<exitCount; iExit++) {
                        cout << table[1].second[iEntrance][iExit] << " ";
                    }
                    cout << endl;
                }
            } else {
                cout << "No other viable connectivity matrices were found." << endl;
            }
        }
    }



    // If the best logP is not good enough, do nothing.
    const double bestLogP = table.front().first;
    if(bestLogP > maxLogP) {
        if(debug) {
            cout << "Not detangling because the best logP is too large." << endl;
        }
        return false;
    }

    // Also check the second best logP.
    if(table.size() >= 2) {
        const double secondBestLogP = table[1].first;
        const double logPDelta = secondBestLogP - bestLogP;
        if(logPDelta < minLogPDelta) {
            if(debug) {
                cout << "Not detangling because the second best logP is too small." << endl;
            }
            return false;
        }

    }

    return false;
}
