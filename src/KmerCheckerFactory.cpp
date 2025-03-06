#include "KmerCheckerFactory.hpp"
#include "Kmer.hpp"
#include "HashedKmerChecker.hpp"
#include "Reads.hpp"
using namespace shasta;



std::shared_ptr<KmerChecker> KmerCheckerFactory::createNew(
    uint64_t k,
    double markerDensity,
    const MappedMemoryOwner& mappedMemoryOwner)
{
    return make_shared<HashedKmerChecker>(
        k,
        markerDensity,
        mappedMemoryOwner);
}



std::shared_ptr<shasta::KmerChecker> KmerCheckerFactory::createFromBinaryData(
    const MappedMemoryOwner& mappedMemoryOwner)
{
    return make_shared<HashedKmerChecker>(mappedMemoryOwner);
}
