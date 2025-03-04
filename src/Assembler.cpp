#include "Assembler.hpp"
#include "AssemblerOptions.hpp"
#include "KmerCheckerFactory.hpp"
#include "MurmurHash2.hpp"
#include "Reads.hpp"
using namespace shasta;

#include "MultithreadedObject.tpp"
template class MultithreadedObject<Assembler>;


// Constructor to be called one to create a new run.
Assembler::Assembler(
    const string& largeDataFileNamePrefixArgument,
    bool createNew,
    size_t largeDataPageSizeArgument) :
    MultithreadedObject(*this)
{
    largeDataFileNamePrefix = largeDataFileNamePrefixArgument;

    if(createNew) {

        // Create a new assembly.
        assemblerInfo.createNew(largeDataName("Info"), largeDataPageSizeArgument);
        assemblerInfo->largeDataPageSize = largeDataPageSizeArgument;
        largeDataPageSize = largeDataPageSizeArgument;

        reads = make_unique<Reads>();
        reads->createNew(
            largeDataName("Reads"),
            largeDataName("ReadNames"),
            largeDataName("ReadMetaData"),
            largeDataName("ReadIdsSortedByName"),
            largeDataPageSize
        );
        // cout << "Created a new assembly with page size " << largeDataPageSize << endl;

    } else {

        // Access an existing assembly.
        assemblerInfo.accessExistingReadWrite(largeDataName("Info"));
        largeDataPageSize = assemblerInfo->largeDataPageSize;

        reads = make_unique<Reads>();
        reads->access(
            largeDataName("Reads"),
            largeDataName("ReadNames"),
            largeDataName("ReadMetaData"),
            largeDataName("ReadIdsSortedByName")
        );
        // cout << "Accessed an existing assembly with page size " << largeDataPageSize << endl;

    }
    SHASTA_ASSERT(largeDataPageSize == assemblerInfo->largeDataPageSize);

    // In both cases, assemblerInfo, reads, readNames, readRepeatCounts are all open for write.

    fillServerFunctionTable();
}









void Assembler::createKmerChecker(
    const KmersOptions& kmersOptions,
    uint64_t threadCount)
{
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    assemblerInfo->k = kmersOptions.k;
    assemblerInfo->kmerGenerationMethod = kmersOptions.generationMethod;

    kmerChecker = KmerCheckerFactory::createNew(
        kmersOptions,
        threadCount,
        getReads(),
        *this);
}



void Assembler::accessKmerChecker()
{
    kmerChecker = KmerCheckerFactory::createFromBinaryData(
        assemblerInfo->k,
        assemblerInfo->kmerGenerationMethod,
        getReads(),
        *this);
}



// Hash a KmerId in such a way that it has the same hash as its reverse
// complement. This is used by alignment method 3 to downsample markers.
uint32_t Assembler::hashKmerId(KmerId kmerId) const
{
    const uint64_t k = assemblerInfo->k;

    // Construct the k-mer and its reverse complement.
    const Kmer kmer(kmerId, k);
    const Kmer kmerRc = kmer.reverseComplement(k);

    // Compute the id of the reverse complement k-mer.
    const KmerId kmerIdRc = KmerId(kmerRc.id(k));

    // Hash the sum of the two KmerIds.
    // This guarantees that we return the same hash
    // for a k-mer and its reverse complement.
    const KmerId sum = kmerId + kmerIdRc;

    return MurmurHash2(&sum, sizeof(sum), 13477);
}
