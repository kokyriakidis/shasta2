
// Shasta.
#include "Assembler.hpp"
#include "AssemblerOptions.hpp"
#include "ConfigurationTable.hpp"
#include "filesystem.hpp"
#include "mode3-Anchor.hpp"
#include "performanceLog.hpp"
#include "Reads.hpp"
#include "Tee.hpp"
#include "timestamp.hpp"
using namespace shasta;

// Boost libraries.
#include <boost/program_options.hpp>
#include  <boost/chrono/process_cpu_clocks.hpp>

//  Linux.
#include <signal.h>
#include <stdlib.h>
#include <unistd.h>


// Standard library.
#include "chrono.hpp"
#include <filesystem>
#include "iostream.hpp"
#include "iterator.hpp"
#include "stdexcept.hpp"



namespace shasta {
    namespace main {

        void main(int argumentCount, const char** arguments);

        void assemble(
            Assembler&,
            const AssemblerOptions&,
            vector<string> inputNames);

        void setupRunDirectory(
            const string& memoryMode,
            const string& memoryBacking,
            size_t& pageSize,
            string& dataDirectory
            );

        void setupHugePages();
        void segmentFaultHandler(int);
        void setupSegmentFaultHandler();

        // Functions that implement --command keywords
        void assemble(const AssemblerOptions&, int argumentCount, const char** arguments);
        void saveBinaryData(const AssemblerOptions&);
        void cleanupBinaryData(const AssemblerOptions&);
        void explore(const AssemblerOptions&);
        void listCommands();

        const std::set<string> commands = {
            "assemble",
            "saveBinaryData",
            "cleanupBinaryData",
            "explore",
            "listCommands"};

    }

    // This is used to duplicate cout output to stdout.log.
    Tee tee;
    ofstream shastaLog;
}



int main(int argumentCount, const char** arguments)
{
    try {
        shasta::main::main(argumentCount, arguments);
    }

    catch(const boost::program_options::error_with_option_name& e) {
        cout << "Invalid option: " << e.what() << endl;
        return 1;
    }

    catch (const runtime_error& e) {
        cout << timestamp << e.what() << endl;
        return 2;
    }

    catch (const std::bad_alloc& e) {
        cout << timestamp << e.what() << endl;
        cout << "Memory allocation failure." << endl;
        cout << "This assembly requires more memory than available." << endl;
        cout << "Rerun on a larger machine." << endl;
        return 2;
    }

    catch (const exception& e) {
        cout << timestamp << e.what() << endl;
        return 3;
    }

    catch (...) {
        cout << timestamp << "Terminated after catching a non-standard exception." << endl;
        return 4;
    }

    return 0;
}



void shasta::main::segmentFaultHandler(int)
{
    char message[] = "\nA segment fault occurred. Please report it by filing an "
        "issue on the Shasta2 repository and attaching the entire log output. "
        "To file an issue, point your browser to https://github.com/paoloshasta/shasta2/issues\n";
    ::write(fileno(stderr), message, sizeof(message));
    ::_exit(1);
}



void shasta::main::setupSegmentFaultHandler()
{
    struct sigaction action;
    ::memset(&action, 0, sizeof(action));
    action.sa_handler = &segmentFaultHandler;
    sigaction(SIGSEGV, &action, 0);
}



void shasta::main::main(int argumentCount, const char** arguments)
{
    setupSegmentFaultHandler();

    // Parse command line options and the configuration file, if one was specified.
    AssemblerOptions assemblerOptions(argumentCount, arguments);

    // Check that we have a valid command.
    auto it = commands.find(assemblerOptions.commandLineOnlyOptions.command);
    if(it ==commands.end()) {
        const string message = "Invalid command " + assemblerOptions.commandLineOnlyOptions.command;
        listCommands();
        throw runtime_error(message);
    }



    // Execute the requested command.
    if(assemblerOptions.commandLineOnlyOptions.command == "assemble") {
        assemble(assemblerOptions, argumentCount, arguments);
        return;
    } else if(assemblerOptions.commandLineOnlyOptions.command == "cleanupBinaryData") {
        cleanupBinaryData(assemblerOptions);
        return;
    } else if(assemblerOptions.commandLineOnlyOptions.command == "saveBinaryData") {
        saveBinaryData(assemblerOptions);
        return;
    } else if(assemblerOptions.commandLineOnlyOptions.command == "explore") {
        explore(assemblerOptions);
        return;
    } else if(assemblerOptions.commandLineOnlyOptions.command == "listCommands") {
        listCommands();
        return;
    }

    // We already checked for a valid command above, so if we get here
    // the above logic is missing code for one of the valid commands.
    SHASTA_ASSERT(0);

}



// Implementation of --command assemble.
void shasta::main::assemble(
    const AssemblerOptions& assemblerOptions,
    int argumentCount, const char** arguments)
{
    SHASTA_ASSERT(assemblerOptions.commandLineOnlyOptions.command == "assemble");


    // Various checks for option validity.

    if(assemblerOptions.commandLineOnlyOptions.configName.empty()) {
        cout <<
            "Option \"--config\" is missing and is now required to "
            "run an assembly.\n"
            "It must specify either a configuration file\n"
            "or one of the following built-in configurations:\n";
        for(const auto& p: configurationTable) {
            cout << p.first << endl;
        }
        throw runtime_error(
            "Option \"--config\" is missing "
            "and is now required to run an assembly.");
    }

    // Check --Kmers.k.
    if(assemblerOptions.kmersOptions.k > 62 or assemblerOptions.kmersOptions.k < 6) {
        throw runtime_error("Invalid value specified for --Kmers.k. "
            "Must be between 6 and 62");
    }
    if((assemblerOptions.kmersOptions.k % 2) == 1) {
        throw runtime_error("Invalid value specified for --Kmers.k. Must be even.");
    }

    // Check that we have at least one input file.
    if(assemblerOptions.commandLineOnlyOptions.inputFileNames.empty()) {
        cout << assemblerOptions.allOptionsDescription << endl;
        throw runtime_error("Specify at least one input file "
            "using command line option \"--input\".");
    }



    // Find absolute paths of the input files.
    // We will use them below after changing directory to the output directory.
    vector<string> inputFileAbsolutePaths;
    for(const string& inputFileName: assemblerOptions.commandLineOnlyOptions.inputFileNames) {
        if(!std::filesystem::exists(inputFileName)) {
            throw runtime_error("Input file not found: " + inputFileName);
        }
        if(!std::filesystem::is_regular_file(inputFileName)) {
            throw runtime_error("Input file is not a regular file: " + inputFileName);
        }
        inputFileAbsolutePaths.push_back(filesystem::getAbsolutePath(inputFileName));
    }



    // Create the assembly directory. If it exists, stop.
    bool exists = std::filesystem::exists(assemblerOptions.commandLineOnlyOptions.assemblyDirectory);
    if (exists) {
        throw runtime_error(
            assemblerOptions.commandLineOnlyOptions.assemblyDirectory +
            " already exists. Remove it first \n"
            "or use --assemblyDirectory to specify a different assembly directory."
        );
    } else {
        SHASTA_ASSERT(std::filesystem::create_directory(assemblerOptions.commandLineOnlyOptions.assemblyDirectory));
    }

    // Make the assembly directory current.
    std::filesystem::current_path(assemblerOptions.commandLineOnlyOptions.assemblyDirectory);

    // Open the performance log.
    openPerformanceLog("performance.log");
    performanceLog << timestamp << "Assembly begins." << endl;

    // Open stdout.log and "tee" (duplicate) stdout to it.
    if(not assemblerOptions.commandLineOnlyOptions.suppressStdoutLog) {
        shastaLog.open("stdout.log");
        tee.duplicate(cout, shastaLog);
    }

    // Echo out the command line options.
    cout << timestamp << "Assembly begins.\nCommand line:" << endl;
    for(int i=0; i<argumentCount; i++) {
        cout << arguments[i] << " ";
    }
    cout << endl;

    // Set up the run directory as required by the memoryMode and memoryBacking options.
    size_t pageSize = 0;
    string dataDirectory;
    setupRunDirectory(
        assemblerOptions.commandLineOnlyOptions.memoryMode,
        assemblerOptions.commandLineOnlyOptions.memoryBacking,
        pageSize,
        dataDirectory);

    // Create the Assembler.
    Assembler assembler(dataDirectory, true, pageSize);

    // Run the assembly.
    assemble(assembler, assemblerOptions, inputFileAbsolutePaths);

    cout << timestamp << "Assembly ends." << endl;
}



// Set up the run directory as required by the memoryMode and memoryBacking options.
void shasta::main::setupRunDirectory(
    const string& memoryMode,
    const string& memoryBacking,
    size_t& pageSize,
    string& dataDirectory
    )
{

    if(memoryMode == "anonymous") {

        if(memoryBacking == "disk") {

            // This combination is meaningless.
            throw runtime_error("\"--memoryMode anonymous\" is not allowed in combination "
                "with \"--memoryBacking disk\".");

        } else if(memoryBacking == "4K") {

            // Anonymous memory on 4KB pages.
            // This combination is the default.
            // It does not require root privilege.
            dataDirectory = "";
            pageSize = 4096;

        } else if(memoryBacking == "2M") {

            // Anonymous memory on 2MB pages.
            // This may require root privilege, which is obtained using sudo
            // and may result in a password prompting depending on sudo set up.
            // Root privilege is not required if 2M pages have already
            // been set up as required.
            setupHugePages();
            pageSize = 2 * 1024 * 1024;

        } else {
            throw runtime_error("Invalid value specified for --memoryBacking: " + memoryBacking +
                "\nValid values are: disk, 4K, 2M.");
        }

    } else if(memoryMode == "filesystem") {

        if(memoryBacking == "disk") {

            // Binary files on disk.
            // This does not require root privilege.
            SHASTA_ASSERT(std::filesystem::create_directory("Data"));
            dataDirectory = "Data/";
            pageSize = 4096;

        } else if(memoryBacking == "4K") {

            // Binary files on the tmpfs filesystem
            // (filesystem in memory backed by 4K pages).
            // This requires root privilege, which is obtained using sudo
            // and may result in a password prompting depending on sudo set up.
            SHASTA_ASSERT(std::filesystem::create_directory("Data"));
            dataDirectory = "Data/";
            pageSize = 4096;
            const string command = "sudo mount -t tmpfs -o size=0 tmpfs Data";
            const int errorCode = ::system(command.c_str());
            if(errorCode != 0) {
                throw runtime_error("Error " + to_string(errorCode) + ": " + strerror(errorCode) +
                    " running command: " + command);
            }

        } else if(memoryBacking == "2M") {

            // Binary files on the hugetlbfs filesystem
            // (filesystem in memory backed by 2M pages).
            // This requires root privilege, which is obtained using sudo
            // and may result in a password prompting depending on sudo set up.
            setupHugePages();
            SHASTA_ASSERT(std::filesystem::create_directory("Data"));
            dataDirectory = "Data/";
            pageSize = 2 * 1024 * 1024;
            const uid_t userId = ::getuid();
            const gid_t groupId = ::getgid();
            const string command = "sudo mount -t hugetlbfs -o pagesize=2M"
                ",uid=" + to_string(userId) +
                ",gid=" + to_string(groupId) +
                " none Data";
            const int errorCode = ::system(command.c_str());
            if(errorCode != 0) {
                throw runtime_error("Error " + to_string(errorCode) + ": " + strerror(errorCode) +
                    " running command: " + command);
            }

        } else {
            throw runtime_error("Invalid value specified for --memoryBacking: " + memoryBacking +
                "\nValid values are: disk, 4K, 2M.");
        }

    } else {
        throw runtime_error("Invalid value specified for --memoryMode: " + memoryMode +
            "\nValid values are: anonymous, filesystem.");
    }
}



// This runs the entire assembly, under the following assumptions:
// - The current directory is the run directory.
// - The Data directory has already been created and set up, if necessary.
// - The input file names are either absolute,
//   or relative to the run directory, which is the current directory.
void shasta::main::assemble(
    Assembler& assembler,
    const AssemblerOptions& assemblerOptions,
    vector<string> inputFileNames)
{

    // Adjust the number of threads, if necessary.
    uint32_t threadCount = assemblerOptions.commandLineOnlyOptions.threadCount;
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }
    cout << "This assembly will use " << threadCount << " threads." << endl;




    // Add reads from the specified input files.
    performanceLog << timestamp << "Begin loading reads from " << inputFileNames.size() << " files." << endl;
    const auto t0 = steady_clock::now();
    for(const string& inputFileName: inputFileNames) {

        assembler.addReads(
            inputFileName,
            assemblerOptions.readsOptions.minReadLength,
            threadCount);
    }

    if(assembler.getReads().readCount() == 0) {
        throw runtime_error("There are no input reads.");
    }
    assembler.computeReadIdsSortedByName();
    assembler.histogramReadLength("ReadLengthHistogram.csv");

    const auto t1 = steady_clock::now();
    performanceLog << timestamp << "Done loading reads from " << inputFileNames.size() << " files." << endl;
    performanceLog << "Read loading took " << seconds(t1-t0) << "s." << endl;

    // Initialize the KmerChecker, which has the information needed
    // to decide if a k-mer is a marker.
    assembler.createKmerChecker(assemblerOptions.kmersOptions, threadCount);

    // Create the markers.
    assembler.createMarkers(threadCount);

    // Create MarkerKmers.
    assembler.createMarkerKmers(threadCount);

    // Create Anchors.
    shared_ptr<mode3::Anchors> anchorsPointer = make_shared<mode3::Anchors>(
        MappedMemoryOwner(assembler),
        assembler.getReads(),
        assembler.assemblerInfo->k,
        assembler.markers(),
        assembler.markerKmers,
        assemblerOptions.assemblyOptions.mode3Options.minAnchorCoverage,
        assemblerOptions.assemblyOptions.mode3Options.maxAnchorCoverage,
        threadCount);

    // Compute Journeys.
    anchorsPointer->computeJourneys(threadCount);

}




// This function sets nr_overcommit_hugepages for 2MB pages
// to a little below total memory.
// If the setting needs to be modified, it acquires
// root privilege via sudo. This may result in the
// user having to enter a password.
void shasta::main::setupHugePages()
{

    // Get the total memory size.
    const uint64_t totalMemoryBytes = sysconf(_SC_PAGESIZE) * sysconf(_SC_PHYS_PAGES);

    // Figure out how much memory we want to allow for 2MB pages.
    const uint64_t MB = 1024 * 1024;
    const uint64_t GB = MB * 1024;
    const uint64_t maximumHugePageMemoryBytes = totalMemoryBytes - 8 * GB;
    const uint64_t maximumHugePageMemoryHugePages = maximumHugePageMemoryBytes / (2 * MB);

    // Check what we have it set to.
    const string fileName = "/sys/kernel/mm/hugepages/hugepages-2048kB/nr_overcommit_hugepages";
    ifstream file(fileName);
    if(!file) {
        throw runtime_error("Error opening " + fileName + " for read.");
    }
    uint64_t currentValue = 0;
    file >> currentValue;
    file.close();

    // If it's set to at least what we want, don't do anything.
    // When this happens, root access is not required.
    if(currentValue >= maximumHugePageMemoryHugePages) {
        return;
    }

    // Use sudo to set.
    const string command =
        "sudo sh -c \"echo " +
        to_string(maximumHugePageMemoryHugePages) +
        " > " + fileName + "\"";
    const int errorCode = ::system(command.c_str());
    if(errorCode != 0) {
        throw runtime_error("Error " + to_string(errorCode) + ": " + strerror(errorCode) +
            " running command: " + command);
    }

}



// Implementation of --command saveBinaryData.
// This copies Data to DataOnDisk.
void shasta::main::saveBinaryData(
    const AssemblerOptions& assemblerOptions)
{
    SHASTA_ASSERT(assemblerOptions.commandLineOnlyOptions.command == "saveBinaryData");

    // Locate the Data directory.
    const string dataDirectory =
        assemblerOptions.commandLineOnlyOptions.assemblyDirectory + "/Data";
    if(!std::filesystem::exists(dataDirectory)) {
        throw runtime_error(dataDirectory + " does not exist, nothing done.");
    }

    // Check that the DataOnDisk directory does not exist.
    const string dataOnDiskDirectory =
        assemblerOptions.commandLineOnlyOptions.assemblyDirectory + "/DataOnDisk";
    if(std::filesystem::exists(dataOnDiskDirectory)) {
        throw runtime_error(dataOnDiskDirectory + " already exists, nothing done.");
    }

    // Copy Data to DataOnDisk.
    const string command = "cp -rp " + dataDirectory + " " + dataOnDiskDirectory;
    const int errorCode = ::system(command.c_str());
    if(errorCode != 0) {
        throw runtime_error("Error " + to_string(errorCode) + ": " + strerror(errorCode) +
            " running command:\n" + command);
    }
    cout << "Binary data successfully saved." << endl;
}



// Implementation of --command cleanupBinaryData.
void shasta::main::cleanupBinaryData(
    const AssemblerOptions& assemblerOptions)
{
    SHASTA_ASSERT(assemblerOptions.commandLineOnlyOptions.command == "cleanupBinaryData");

    // Locate the Data directory.
    const string dataDirectory =
        assemblerOptions.commandLineOnlyOptions.assemblyDirectory + "/Data";
    if(!std::filesystem::exists(dataDirectory)) {
        cout << dataDirectory << " does not exist, nothing done." << endl;
        return;
    }

    // Unmount it and remove it.
    ::system(("sudo umount " + dataDirectory).c_str());
    const int errorCode = ::system(string("rm -rf " + dataDirectory).c_str());
    if(errorCode != 0) {
        throw runtime_error("Error " + to_string(errorCode) + ": " + strerror(errorCode) +
            " removing " + dataDirectory);
    }
    cout << "Cleanup of " << dataDirectory << " successful." << endl;

    // If the DataOnDisk directory exists, create a symbolic link
    // Data->DataOnDisk.
    const string dataOnDiskDirectory =
        assemblerOptions.commandLineOnlyOptions.assemblyDirectory + "/DataOnDisk";
    if(std::filesystem::exists(dataOnDiskDirectory)) {
        std::filesystem::current_path(assemblerOptions.commandLineOnlyOptions.assemblyDirectory);
        const string command = "ln -s DataOnDisk Data";
        ::system(command.c_str());
    }

}

// Implementation of --command explore.
void shasta::main::explore(
    const AssemblerOptions& assemblerOptions)
{

    // Go to the assembly directory.
    std::filesystem::current_path(assemblerOptions.commandLineOnlyOptions.assemblyDirectory);

    // Check that we have the binary data.
    if(!std::filesystem::exists("Data")) {
        throw runtime_error("Binary directory \"Data\" not available "
        " in assembly directory " +
        assemblerOptions.commandLineOnlyOptions.assemblyDirectory +
        ". Use \"--memoryMode filesystem\", possibly followed by "
        "\"--command saveBinaryData\" and \"--command cleanupBinaryData\" "
        "if you want to make sure the binary data are persistently available on disk."
        );
        return;
    }

    // Create the Assembler.
    Assembler assembler("Data/", false, 0);

    // Access all available binary data.
    assembler.httpServerData.assemblerOptions = &assemblerOptions;
    assembler.accessAllSoft();


    // Start the http server.
    bool localOnly;
    bool sameUserOnly;
    if(assemblerOptions.commandLineOnlyOptions.exploreAccess == "user") {
        localOnly = true;
        sameUserOnly = true;
    } else if(assemblerOptions.commandLineOnlyOptions.exploreAccess == "local") {
        localOnly = true;
        sameUserOnly = false;
    } else if (assemblerOptions.commandLineOnlyOptions.exploreAccess == "unrestricted"){
        localOnly = false;
        sameUserOnly = false;
    } else {
        throw runtime_error("Invalid value specified for --exploreAccess. "
            "Only use this option if you understand its security implications."
        );
    }
    assembler.explore(
        assemblerOptions.commandLineOnlyOptions.port,
        localOnly,
        sameUserOnly);
}



void shasta::main::listCommands()
{
    cout << "Valid commands are:" << endl;
    for(const string& command: commands) {
        cout << command << endl;
    }
}
