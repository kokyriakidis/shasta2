// Shasta.
#include "PermutationDetangler.hpp"
#include "orderPairs.hpp"
#include "Tangle.hpp"
using namespace shasta;

// Boost libraries.
#include <boost/math/distributions/chi_squared.hpp>

// Standard library.
#include "iterator.hpp"
#include "utility.hpp"



PermutationDetangler::PermutationDetangler(
    uint64_t minCommonCoverage,
    const double epsilon,
    const double maxLogP,
    const double minLogPDelta) :
    minCommonCoverage(minCommonCoverage),
    epsilon(epsilon),
    maxLogP(maxLogP),
    minLogPDelta(minLogPDelta)
{}



bool PermutationDetangler::operator(
    )(
    Tangle& tangle)
{

    const TangleMatrix& tangleMatrix = tangle.tangleMatrix;
    const vector< vector<uint64_t> >& matrix = tangleMatrix.tangleMatrix;

    // The PermutationDetangler only works on Tangles where the
    // numbers of Entrances and Exits are the same and at least 2.
    const uint64_t n = tangleMatrix.entrances.size();
    if(n < 2) {
        return false;
    }
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

    // Extract common coverage for each entrance and exit and total common coerage.
    vector<uint64_t> entranceCommonCoverage(n);
    for(uint64_t iEntrance=0; iEntrance<n; iEntrance++) {
        entranceCommonCoverage[iEntrance] = tangleMatrix.entrances[iEntrance].commonCoverage;
    }
    vector<uint64_t> exitCommonCoverage(n);
    for(uint64_t iExit=0; iExit<n; iExit++) {
        exitCommonCoverage[iExit] = tangleMatrix.exits[iExit].commonCoverage;
    }
    const uint64_t totalCommonCoverage = accumulate(entranceCommonCoverage.begin(), entranceCommonCoverage.end(), 0UL);
    SHASTA_ASSERT(totalCommonCoverage == accumulate(exitCommonCoverage.begin(), exitCommonCoverage.end(), 0UL));

    if(false) {
        cout << "Entrance common coverage:";
        for(uint64_t iEntrance=0; iEntrance<n; iEntrance++) {
            cout << " " << entranceCommonCoverage[iEntrance];
        }
        cout << endl;
        cout << "Exit common coverage:";
        for(uint64_t iExit=0; iExit<n; iExit++) {
            cout << " " << exitCommonCoverage[iExit];
        }
        cout << endl;
    }

    // Chi squared distribution used below to test each permutation.
    boost::math::chi_squared_distribution chi2Distribution(double(n * n - 1));

    // Permutations and their logP.
    // Only store the ones for which logP is not inf.
    vector< pair<vector<uint64_t>, double> > permutationTable;

    // Try all possible permutations of the Exits, leaving the Entrances in the same order.
    // For each permutation, compute the permuted tangle matrix and perform a chi square test on it.
    vector<uint64_t> permutation(n);
    std::iota(permutation.begin(), permutation.end(), 0);
    vector< vector<uint64_t> > permutedMatrix(n, vector<uint64_t>(n));
    vector<uint64_t> permutedExitCommonCoverage(n);
    do {
        if(false) {
            cout << "Trying permutation ";
            copy(permutation.begin(), permutation.end(), ostream_iterator<uint64_t>(cout, " "));
            cout << endl;
        }

        // Compute the permuted tangle matrix.
        for(uint64_t iEntrance=0; iEntrance<n; iEntrance++) {
            for(uint64_t iExit=0; iExit<n; iExit++) {
                permutedMatrix[iEntrance][iExit] = matrix[iEntrance][permutation[iExit]];
            }
        }

        if(false) {
            for(uint64_t iEntrance=0; iEntrance<n; iEntrance++) {
                for(uint64_t iExit=0; iExit<n; iExit++) {
                    cout << iEntrance << " " << iExit << " " << permutedMatrix[iEntrance][iExit] << endl;
                }
            }
        }

        // Compute permuted exit common coverage.
        for(uint64_t iExit=0; iExit<n; iExit++) {
            permutedExitCommonCoverage[iExit] = exitCommonCoverage[permutation[iExit]];
        }

        if(false) {
            cout << "Permuted exit common coverage:";
            for(uint64_t iExit=0; iExit<n; iExit++) {
                cout << " " << permutedExitCommonCoverage[iExit];
            }
            cout << endl;
        }



        // Assuming this permutation, create expected values for the tangle matrix
        // under various assumptions.

        // Random.
        vector< vector<double> > tangleMatrixRandom(n, vector<double>(n, 0.));
        for(uint64_t iEntrance=0; iEntrance<n; iEntrance++) {
            for(uint64_t iExit=0; iExit<n; iExit++) {
                tangleMatrixRandom[iEntrance][iExit] =
                    double(entranceCommonCoverage[iEntrance]) *
                    double(permutedExitCommonCoverage[iExit]) /
                    double(totalCommonCoverage);
            }
        }

        // Assuming this permutation, the ideal tangle matrix is diagonal.
        vector< vector<double> > tangleMatrixIdeal(n, vector<double>(n, 0.));
        for(uint64_t i=0; i<n; i++) {
            tangleMatrixIdeal[i][i] = 0.5 * (double(entranceCommonCoverage[i]) + double(permutedExitCommonCoverage[i]));
        }

        // Compute the expected tangle matrix under non-ideal assumptions.
        vector< vector<double> > tangleMatrixNonIdeal(n, vector<double>(n, 0.));
        for(uint64_t iEntrance=0; iEntrance<n; iEntrance++) {
            for(uint64_t iExit=0; iExit<n; iExit++) {
                tangleMatrixNonIdeal[iEntrance][iExit] =
                    epsilon * tangleMatrixRandom[iEntrance][iExit] +
                    (1. - epsilon) * tangleMatrixIdeal[iEntrance][iExit];
            }
        }

        if(false) {
            cout << "Random tangle matrix for this permutation:" << endl;
            for(uint64_t iEntrance=0; iEntrance<n; iEntrance++) {
                for(uint64_t iExit=0; iExit<n; iExit++) {
                    cout << iEntrance << " " << iExit << " " << tangleMatrixRandom[iEntrance][iExit] << endl;
                }
            }

            cout << "Expected ideal tangle matrix for this permutation:" << endl;
            for(uint64_t iEntrance=0; iEntrance<n; iEntrance++) {
                for(uint64_t iExit=0; iExit<n; iExit++) {
                    cout << iEntrance << " " << iExit << " " << tangleMatrixIdeal[iEntrance][iExit] << endl;
                }
            }

            cout << "Expected non-ideal tangle matrix for this permutation:" << endl;
            for(uint64_t iEntrance=0; iEntrance<n; iEntrance++) {
                for(uint64_t iExit=0; iExit<n; iExit++) {
                    cout << iEntrance << " " << iExit << " " << tangleMatrixNonIdeal[iEntrance][iExit] << endl;
                }
            }
        }

        // Do a chi-square test.
        double chi2 = 0.;
        for(uint64_t iEntrance=0; iEntrance<n; iEntrance++) {
            for(uint64_t iExit=0; iExit<n; iExit++) {
                const double expected = tangleMatrixNonIdeal[iEntrance][iExit];
                const double delta = double(permutedMatrix[iEntrance][iExit]) - expected;
                chi2 += delta * delta / expected;
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

            // Store this permutation and its logP.
            permutationTable.push_back(make_pair(permutation, logP));

            if(false) {
                cout << "Chi square for this permutation: " << chi2 <<
                    ", logP " << logP << " dB." << endl;
            }
        }


    } while(std::next_permutation(permutation.begin(), permutation.end()));

    // If we did not find any usable permutation, do nothing.
    if(permutationTable.empty()) {
        if(debug) {
            cout << "The permutation table is empty." << endl;
        }
        return false;
    }

    // Sort the permutation table by increasing logP.
    sort(permutationTable.begin(), permutationTable.end(),
        OrderPairsBySecondOnly<vector<uint64_t>, double>());

    if(debug) {
        cout << "Permutation table:" << endl;
        for(const auto& p: permutationTable) {
            const vector<uint64_t>& permutation = p.first;
            const double logP = p.second;

            for(uint64_t i=0; i<n; i++) {
                cout << "(" <<
                    i << "," <<
                    permutation[i] << ") ";
            }
            cout << logP << endl;
        }
    }

    // If the best logP is not good enough, do nothing.
    const double bestLogP = permutationTable.front().second;
    if(bestLogP > maxLogP) {
        if(debug) {
            cout << "Not detangling because the best logP is too large." << endl;
        }
        return false;
    }

    // Also check the second best logP.
    if(permutationTable.size() >= 2) {
        const double secondBestLogP = permutationTable[1].second;
        const double logPDelta = secondBestLogP - bestLogP;
        if(logPDelta < minLogPDelta) {
            if(debug) {
                cout << "Not detangling because the second best logP is too small." << endl;
            }
            return false;
        }

    }

    if(debug) {
        cout << "This superbubble would be detangled." << endl;
    }


    // Detangle using the best permutation.
    const vector<uint64_t>& bestPermutation = permutationTable.front().first;
    for(uint64_t i=0; i<n; i++) {
        tangle.connect(i, bestPermutation[i]);
    }
    tangle.detangle();

    return true;

}
