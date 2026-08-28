#pragma once

// Shasta
#include "invalid.hpp"

// Standard library.
#include "cstdint.hpp"
#include "iosfwd.hpp"
#include "vector.hpp"



// Likelihood ratio test of the tangle matrix (G test).
// https://en.wikipedia.org/wiki/G-test
namespace shasta2 {
    class GTest;
}



class shasta2::GTest {
public:
    GTest(
        const vector< vector<uint64_t> >& tangleMatrix,
        double epsilon,
        bool onlyConsiderInjective,
        bool onlyConsiderPermutation);
    GTest(
        const vector< vector<double> >& tangleMatrix,
        double epsilon,
        bool onlyConsiderInjective,
        bool onlyConsiderPermutation);
    bool success = false;

    void writeHtml(ostream&) const;

    class Hypothesis {
    public:
        vector< vector<bool> > connectivityMatrix;
        double G = invalid<double>;

        Hypothesis(
            const vector< vector<bool> >& connectivityMatrix,
            double G) :
            connectivityMatrix(connectivityMatrix),
            G(G)
            {}

        Hypothesis() {}


        // Sort by G.
        bool operator<(const Hypothesis& that) const {
            return G < that.G;
        }
    };
    vector<Hypothesis> hypotheses;

    // Return true if there is a single exit for each entrance.
    static bool isForwardInjective(const vector< vector<bool> >& connectivityMatrix);

    // Return true if there is a single entrance for exit entrance.
    static bool isBackwardInjective(const vector< vector<bool> >& connectivityMatrix);

    // The GTest is considered positive if the following is true:
    // - The GTest ran successfully (that is, success is set to true).
    // - The G value of the first hypothesis is no more than maxLogP.
    // - If a second hypothesis is present, the G difference between
    //   the second hypothesis and the first is at least minLogP.
    // In this case, the hypotheses vector is guaranteed to not be empty,
    // and the first hypothesis in the vector is likely to be reliable.
    bool isPositive(
        double maxLogP,
        double minLogPDelta) const;

    // Summarize the GTest results.
    void summarize(ostream&) const;

private:
    void run(
        const vector< vector<double> >& tangleMatrix,
        double epsilon,
        bool onlyConsiderInjective,
        bool onlyConsiderPermutation);

};
