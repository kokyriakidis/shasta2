// Shasta.
#include "PermutationDetangler.hpp"
#include "orderPairs.hpp"
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



bool PermutationDetangler::operator()(Tangle&)
{
    return false;
}
