#include "KmerCheckerFactory.hpp"
#include "Kmer.hpp"
#include "HashedKmerChecker.hpp"
#include "Reads.hpp"
using namespace shasta2;



std::shared_ptr<KmerChecker> KmerCheckerFactory::createNew(
    uint64_t k,
    double markerDensity)
{
    return make_shared<HashedKmerChecker>(
        k,
        markerDensity);
}
