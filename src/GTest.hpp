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

private:
    void run(
        const vector< vector<double> >& tangleMatrix,
        double epsilon,
        bool onlyConsiderInjective,
        bool onlyConsiderPermutation);

};
