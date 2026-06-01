#pragma once
#include "LongBaseSequence.hpp"
#include "MultithreadedObject.hpp"

#include "string.hpp"


namespace shasta2 {
    class FastaLoader;

    class Reads;
}



class shasta2::FastaLoader : public MultithreadedObject<FastaLoader> {
public:

    FastaLoader(
        uint64_t minReadLength,
        size_t threadCount,
        Reads& reads);

    ~FastaLoader();

    void read(const string& fileName);

private:

    uint64_t minReadLength;
    size_t threadCount;
    Reads& reads;

    FastaLoader(const FastaLoader&) = delete;
    FastaLoader& operator=(const FastaLoader&) const = delete;

    // A page-aligned buffer.
    const uint64_t bufferSize = 256 * 1024 * 1024;
    char* buffer = 0;
    void allocateBuffer();
    void freeBuffer();

    int fileDescriptor = 0;
    bool endOfFile = false;
    void openFile(const string& fileName);
    void closeFile();

    void readBlock();

    // The size of the last block we read.
    // It will always be equal to bufferSize, except when we
    // reach the end of the file, in which case it could be less.
    uint64_t blockSize = 0;

    // The file offset of the current block.
    uint64_t blockOffset = 0;

    bool lastCharacterOfPreviousBlockWasEndLine = false;

    // All the places in the current block where a new read begins.
    // This happens at '>' characters that are preceded by '\n',
    // plus the mandatory '>' at the beginning of the file.
    vector<uint64_t> readBeginPositions;
    void findReadBeginPositions();


    // Each thread started by runThreads stores the reads it finds.
    // These are reads entirely contained in a block.
    // The main thread stores the remaining reads.
    class ReadInfo {
    public:

        // The file offset where this read begins.
        uint64_t fileOffset;

        string name;
        LongBaseSequence sequence;

        ReadInfo(
            const std::string_view&,
            uint64_t fileOffset,
            vector<Base>&);
    };

    // This has size threadCount+1. The position at index threadCount
    // is used by the main thread.
    vector< vector<ReadInfo> > readInfos;

    void threadFunction(uint64_t threadId);

};
