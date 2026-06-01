// Shasta2.
#include "FastaLoader.hpp"
#include "invalid.hpp"
#include "Reads.hpp"
using namespace shasta2;

// Linux.
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

// Explicit instantiation.
#include "MultithreadedObject.tpp"
template class MultithreadedObject<FastaLoader>;



FastaLoader::FastaLoader(
    uint64_t minReadLength,
    size_t threadCount,
    Reads& reads) :
    MultithreadedObject(*this),
    minReadLength(minReadLength),
    threadCount(threadCount),
    reads(reads)
{

}



FastaLoader::~FastaLoader()
{
    freeBuffer();
    closeFile();
}



void FastaLoader::read(const string& fileName)
{
    const bool debug = false;

    allocateBuffer();
    openFile(fileName);

    readInfos.clear();
    readInfos.resize(threadCount + 1);

    // A string to contain the partial read in the previous block.
    string partialRead;
    uint64_t partialReadFileOffset = invalid<uint64_t>;

    vector<Base> sequence;



    // Main loop. Each iteration reads and processes one block.
    blockOffset = 0;
    while(true) {

        // Read a block.
        readBlock();
        if(blockSize == 0) {
            break;
        }
        if(debug) {
            cout << timestamp << "Read " << blockSize << " bytes from fasta file." << endl;
            cout << "    This block covers file offsets " << blockOffset << " to " <<
                blockOffset + blockSize << endl;
        }

        // Find all the places in this block where a new read begins.
        // This happens at '>' characters that are preceded by '\n',
        // plus the mandatory '>' at the beginning of the file.
        findReadBeginPositions();
        if(debug) {
            cout << readBeginPositions.size() << " reads begin in this block." << endl;
        }

        // Use multithreaded code to parse the reads entirely contained in this block.
        // Each thread stores in its own vector<ReadInfo>.
        if(readBeginPositions.size() > 1) {
            const uint64_t batchSize = 10;
            setupLoadBalancing(readBeginPositions.size() - 1, batchSize);
            runThreads(&FastaLoader::threadFunction, threadCount);
        }

        // Special case when no read begins in this block.
        // This should almost never happen.
        if(readBeginPositions.empty()) {
            partialRead.append(buffer, blockSize);
        } else {

            // At least one read begins in this block.
            if(blockOffset > 0) {
                if(partialReadFileOffset == invalid<uint64_t>) {
                    throw runtime_error("Malformed fasta file.");
                }
                const uint64_t initialPartialReadBegin = 0;
                const uint64_t initialPartialReadEnd = readBeginPositions.front();
                partialRead.append(buffer + initialPartialReadBegin, initialPartialReadEnd - initialPartialReadBegin);
                readInfos.back().emplace_back(partialRead, partialReadFileOffset, sequence);
                if(readInfos.back().back().sequence.baseCount < minReadLength) {
                    readInfos.back().pop_back();
                }
            }

            // Handle the partial read at the end of this clock.
            const uint64_t finalPartialReadBegin = readBeginPositions.back();
            const uint64_t finalPartialReadEnd = blockSize;
            partialRead = string(buffer + finalPartialReadBegin, finalPartialReadEnd - finalPartialReadBegin);
            partialReadFileOffset = blockOffset + finalPartialReadBegin;
        }


        if(endOfFile) {
            break;
        }

        // Store the file offset for the next block.
        blockOffset += blockSize;
        lastCharacterOfPreviousBlockWasEndLine = (buffer[blockSize-1] == '\n');
    }

    // Handle the partial read at the end of the file.
    readInfos.back().emplace_back(partialRead, partialReadFileOffset, sequence);
    if(readInfos.back().back().sequence.baseCount < minReadLength) {
        readInfos.back().pop_back();
    }


    // Now we have to copy the reads from the ReadInfos to the Reads.
    // We want to store them in the same order as they appear in the file,
    // so we have to sort them by fileOffset.
    class Info {
    public:
        uint64_t threadId;
        uint64_t indexInThread;
        uint64_t fileOffset;
        bool operator<(const Info& that) const {
            return fileOffset < that.fileOffset;
        }
    };
    vector<Info> table;
    for(uint64_t threadId=0; threadId<=threadCount; threadId++) {
        const vector<ReadInfo>& threadReadInfos = readInfos[threadId];
        for(uint64_t indexInThread=0; indexInThread<threadReadInfos.size(); indexInThread++) {
            const ReadInfo& readInfo = threadReadInfos[indexInThread];
            table.emplace_back(Info({threadId, indexInThread, readInfo.fileOffset}));
        }
    }
    sort(table.begin(), table.end());
    for(const Info& info: table) {
        const ReadInfo& readInfo = readInfos[info.threadId][info.indexInThread];
        reads.readNames.appendVector(readInfo.name.begin(), readInfo.name.end());
        reads.reads.append(readInfo.sequence);
    }


    readInfos.clear();

}



// Find all the places in the current block where a new read begins.
// This happens at '>' characters that are preceded by '\n',
// plus the mandatory '>' at the beginning of the file.
void FastaLoader::findReadBeginPositions()
{
    readBeginPositions.clear();

    if(blockOffset == 0) {
        // The first block must begin with '>'.
        if(buffer[0] != '>') {
            throw runtime_error("Fasta file does not begin with \">\".");
        }
        readBeginPositions.push_back(0);
    } else {
        if(lastCharacterOfPreviousBlockWasEndLine and (buffer[0] == '>')) {
            readBeginPositions.push_back(0);
        }
    }
    const std::string_view s(buffer, blockSize);
    uint64_t position = 1;
    while(true) {
        position = s.find_first_of('>', position);
        if(position == string::npos) {
            break;
        }
        SHASTA2_ASSERT(buffer[position] == '>');
        if(buffer[position-1] == '\n') {
            readBeginPositions.push_back(position);
        }
        ++position;
    }
}



void FastaLoader::openFile(const string& fileName)
{
    closeFile();

    // Open the file, using O_DIRECT if possible.
    // This is done so the fasta file does not go in the Linux cache,
    // eliminating the overhead to later remove the pages
    // from the cache when needed by the assembly.
    // This is especially useful when using large pages
    // (--memory-mode filesystem --memory-backing 2M).
    fileDescriptor = ::open(fileName.c_str(), O_RDONLY | O_DIRECT);
    if(fileDescriptor <= 0) {
        // O_DIRECT failed, try without it.
        fileDescriptor = ::open(fileName.c_str(), O_RDONLY);
    }
    if(fileDescriptor <= 0) {
        throw runtime_error("Error reading " + fileName);
    }
    endOfFile = false;
    blockOffset = 0;
    lastCharacterOfPreviousBlockWasEndLine = false;
}



void FastaLoader::allocateBuffer()
{
    if(not buffer) {
        void* pointer = ::mmap(
            0,
            bufferSize,
            PROT_READ | PROT_WRITE,
            MAP_PRIVATE | MAP_ANONYMOUS,
            -1,
            0);
        if(pointer == MAP_FAILED) {
            throw runtime_error("Error allocating FastaLoader buffer.");
        }
        buffer = static_cast<char*>(pointer);
    }

}



void FastaLoader::freeBuffer()
{
    if(buffer) {
        ::munmap(buffer, bufferSize);
    }
}



void FastaLoader::closeFile()
{
    if(fileDescriptor > 0) {
        ::close(fileDescriptor);
    }
}



void FastaLoader::readBlock()
{
    // Try to read bufferSize bytes from the currently open file and
    // store in blockSize the number of bytes read.
    // That will always be equal to bufferSize except at
    // the end of file.

    uint64_t bytesToRead = bufferSize;
    char* bufferPointer = buffer;
    blockSize = 0;
    while(bytesToRead) {
        const int64_t bytesRead = ::read(fileDescriptor, bufferPointer, bytesToRead);
        if(bytesRead == -1) {
            throw runtime_error("Error reading fasta file.");
        }
        if(bytesRead == 0) {
            // End of file.
            endOfFile = true;
            return;
        }
        bufferPointer += bytesRead;
        blockSize += bytesRead;
        bytesToRead -= bytesRead;
    }
}



void FastaLoader::threadFunction(uint64_t threadId)
{

    // This thread will store here the reads it finds.
    vector<ReadInfo>& treadReadInfo = readInfos[threadId];

    vector<Base> sequence;

    // Loop over all batches assigned to this thread.
    uint64_t begin;
    uint64_t end;
    while(getNextBatch(begin, end)) {

        // Loop over the reads in this batch.
        for(uint64_t i=begin; i!=end; ++i) {


            // Create a string_view with the characters
            // that define this read.
            char const* readBegin = buffer + readBeginPositions[i];
            char const* readEnd = buffer + readBeginPositions[i+1];
            const std::string_view s(readBegin, readEnd);

            // Use the string_view to create a new ReadInfo and store it.
            treadReadInfo.emplace_back(s, blockOffset + readBeginPositions[i], sequence);
            if(treadReadInfo.back().sequence.baseCount < minReadLength) {
                treadReadInfo.pop_back();
            }
        }

    }
}



FastaLoader::ReadInfo::ReadInfo(
    const std::string_view& s,
    uint64_t fileOffset,
    vector<Base>& sequenceVector) :
    fileOffset(fileOffset)
{
    // Because of how we constructed readBeginPositions,
    // the first character is guaranteed t be '>'
    // and the last character is guaranteed to be '\n'.
    SHASTA2_ASSERT(s.front() == '>');
    SHASTA2_ASSERT(s.back() == '\n');

    // Find the line-end which marks the end of the header line.
    const uint64_t lineEndPosition = s.find_first_of('\n');

    // Create the header, without the initial '>'.
    const std::string_view header(s.begin()+1, s.begin() + lineEndPosition);

    // Store the name.
    const uint64_t firstBlankPosition = s.find_first_of(' ');
    name = header.substr(1, firstBlankPosition - 1);

    // Store the sequence.
    sequenceVector.clear();
    const uint64_t sequenceBegin = lineEndPosition + 1;
    const uint64_t sequenceEnd = s.size();
    for(uint64_t i=sequenceBegin; i!=sequenceEnd; i++) {
        const char c = s[i];
        const Base b = Base::fromCharacterNoException(c);

        if(b.isValid()) {
            sequenceVector.push_back(b);
        } else {
            if((c != ' ') and (c != '\n')) {
                throw runtime_error("Invalid base character " + string(1, c));
            }
        }
    }
    sequence.assign(sequenceVector);
    sequenceVector.clear();

}
