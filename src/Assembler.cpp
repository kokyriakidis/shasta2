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
            largeDataName("ReadIdsSortedByName")
        );

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

    kmerChecker = KmerCheckerFactory::createNew(
        kmersOptions,
        threadCount,
        getReads(),
        *this);
}



void Assembler::accessKmerChecker()
{
    kmerChecker = KmerCheckerFactory::createFromBinaryData(*this);
}
