#include "KmerCheckerFactory.hpp"
#include "Kmer.hpp"
#include "HashedKmerChecker.hpp"
#include "AssemblerOptions.hpp"
#include "Reads.hpp"
using namespace shasta;



std::shared_ptr<KmerChecker> KmerCheckerFactory::createNew(
    const KmersOptions& kmersOptions,
    uint64_t /* threadCount */,
    const Reads& /* reads */,
    const MappedMemoryOwner& mappedMemoryOwner)
{
    return make_shared<HashedKmerChecker>(
        kmersOptions.k,
        kmersOptions.probability,
        mappedMemoryOwner);
}



std::shared_ptr<shasta::KmerChecker> KmerCheckerFactory::createFromBinaryData(
    const MappedMemoryOwner& mappedMemoryOwner)
{
    return make_shared<HashedKmerChecker>(mappedMemoryOwner);
}
