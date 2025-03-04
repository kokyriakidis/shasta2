#include "AssemblerOptions.hpp"
#include "ConfigurationTable.hpp"
#include "filesystem.hpp"
using namespace shasta;

// Boost libraries.
#include <boost/tokenizer.hpp>

// Standard library.
#include "fstream.hpp"
#include "stdexcept.hpp"



// Constructor from a command line.
// If the command line includes a --config option,
// the specified built-in configuration or configuration file
// is used to fill the AssemblyOptions,
// but values specified on the command line take precedence.
AssemblerOptions::AssemblerOptions(int argumentCount, const char** arguments) :
    commandLineOnlyOptionsDescription("Options allowed only on the command line"),
    configurableOptionsDescription("Options allowed on the command line and in the config file")
{
    using boost::program_options::positional_options_description;
    using boost::program_options::value;
    using boost::program_options::variables_map;
    using boost::program_options::command_line_parser;

    addCommandLineOnlyOptions();
    addConfigurableOptions();

    // Define allOptionsDescription as the union of
    // commandLineOnlyOptionsDescription and configurableOptionsDescription.
    allOptionsDescription.add(commandLineOnlyOptionsDescription);
    allOptionsDescription.add(configurableOptionsDescription);

    // Machinery to capture invalid positional options.
    allOptionsIncludingInvalidDescription.add(allOptionsDescription);
    allOptionsIncludingInvalidDescription.add_options()
        ("invalidOption",
            value< vector<string> >(&invalidPositionalOptions));
    positional_options_description positionalOptionsDescription;
    positionalOptionsDescription.add("invalidOption", -1);


    // Get options from the command line.
    // These take precedence over values entered in the config file.
    variables_map variablesMap;
    store(command_line_parser(argumentCount, arguments).
          options(allOptionsIncludingInvalidDescription).
          positional(positionalOptionsDescription).
          run(),
          variablesMap);
    notify(variablesMap);

    // If any invalid positional options were specified,
    // stop here.
    if(!invalidPositionalOptions.empty()) {
        string message = "Positional options are not allowed. "
            "The following positional options were used:";
        for(const string& s: invalidPositionalOptions) {
            message += " ";
            message += s;
        }
        throw runtime_error(message);
    }

    // If help was requested, write the help message and exit.
    if (variablesMap.count("help")) {
        cout << allOptionsDescription << endl;
        ::exit(0);
    }

    // If the command is explore and the configuration is mnot specified, use the
    // configuration file in the assembly directory.
    if((commandLineOnlyOptions.command == "explore") and commandLineOnlyOptions.configName.empty()) {
        commandLineOnlyOptions.configName = commandLineOnlyOptions.assemblyDirectory + "/shasta.conf";
    }


    // Get options from the config file or built-in configuration, if one was specified.
    if(!commandLineOnlyOptions.configName.empty()) {

        // Look it up in our configuration table.
        const string* builtInConfiguration = getConfiguration(commandLineOnlyOptions.configName);

        if(builtInConfiguration) {

            // --config specifies a built-in configuration. Get options from that.
            std::istringstream s(*builtInConfiguration);
            store(parse_config_file(s, configurableOptionsDescription), variablesMap);
            notify(variablesMap);

        } else {

            // See if --config specifies a configuration file.
            ifstream configFile(commandLineOnlyOptions.configName);
            if (!configFile) {
                string message =
                    "The --config option does not specify a "
                    "valid built-in configuration or configuration file. "
                    "Valid built-in configurations are the following:\n";
                for(const auto& p: configurationTable) {
                    message += p.first + "\n";
                }
                throw runtime_error(message);
            }
            store(parse_config_file(configFile, configurableOptionsDescription), variablesMap);
            notify(variablesMap);
        }
    }
}



// Constructor from a configuration file.
// This only fills in the configurable options specified in
// the given configuration file. Command line only options
// are left at their defaults.
AssemblerOptions::AssemblerOptions(const string& fileName)
{

    using boost::program_options::positional_options_description;
    using boost::program_options::value;
    using boost::program_options::variables_map;
    using boost::program_options::command_line_parser;

    addConfigurableOptions();

    ifstream configFile(fileName);
    if (not configFile) {
        throw runtime_error("Invalid configuration file " + fileName + " specified.\n");
    }
    variables_map variablesMap;
    store(parse_config_file(configFile, configurableOptionsDescription), variablesMap);
    notify(variablesMap);

}



// Add non-configurable options to the Boost option description object.
// These are options that can only be used on the command line
// (not in the configuration file).
void AssemblerOptions::addCommandLineOnlyOptions()
{
    using boost::program_options::value;
    using boost::program_options::bool_switch;

    commandLineOnlyOptionsDescription.add_options()

        ("help,h",
        "Write a help message.")

        ("version,v",
        "Identify the Shasta version.")

        ("config",
        value<string>(&commandLineOnlyOptions.configName),
        "Configuration name. Can be the name of a built-in configuration "
        "or the name of a configuration file.")

        ("input",
        value< vector<string> >(&commandLineOnlyOptions.inputFileNames)->multitoken(),
        "Names of input files containing reads. Specify at least one.")

        ("anchors",
        value< vector<string> >(&commandLineOnlyOptions.anchorFileNames)->multitoken(),
        "Names of input files containing anchors for mode 3 assembly.")

        ("assemblyDirectory",
        value<string>(&commandLineOnlyOptions.assemblyDirectory)->
        default_value("ShastaRun"),
        "Name of the output directory. If command is assemble, this directory must not exist.")

        ("command",
        value<string>(&commandLineOnlyOptions.command)->
        default_value("assemble"),
        "Command to run. Must be one of: "
        "assemble, saveBinaryData, cleanupBinaryData, explore, createBashCompletionScript")

        ("memoryMode",
        value<string>(&commandLineOnlyOptions.memoryMode)->
        default_value("anonymous"),
        "Specify whether allocated memory is anonymous or backed by a filesystem. "
        "Allowed values: anonymous, filesystem.")

        ("memoryBacking",
        value<string>(&commandLineOnlyOptions.memoryBacking)->
        default_value("4K"),
        "Specify the type of pages used to back memory.\n"
        "Allowed values: disk, 4K , 2M (for best performance). "
        "All combinations (memoryMode, memoryBacking) are allowed "
        "except for (anonymous, disk).\n"
        "Some combinations require root privilege, which is obtained using sudo "
        "and may result in a password prompting depending on your sudo set up.")

        ("threads",
        value<uint32_t>(&commandLineOnlyOptions.threadCount)->
        default_value(0),
        "Number of threads, or 0 to use one thread per virtual processor.")
        
        ("suppressStdoutLog",
        bool_switch(&commandLineOnlyOptions.suppressStdoutLog)->
        default_value(false),
        "Suppress echoing stdout to stdout.log.")

        ("exploreAccess",
        value<string>(&commandLineOnlyOptions.exploreAccess)->
        default_value("user"),
        "Specify allowed access for --command explore. "
        "Allowed values: user, local, unrestricted. "
        "DO NOT CHANGE FROM DEFAULT VALUE WITHOUT UNDERSTANDING THE "
        "SECURITY IMPLICATIONS."
        )

        ("port",
        value<uint16_t>(&commandLineOnlyOptions.port)->
        default_value(17100),
        "Port to be used by the http server (command --explore).")

        ("alignmentsPafFile",
        value<string>(&commandLineOnlyOptions.alignmentsPafFile),
        "The name of a PAF file containing alignments of reads to "
        "a reference. Only used for --command explore, for display of the alignment "
        "candidate graph. Experimental."
        )

        ("saveBinaryData",
        bool_switch(&commandLineOnlyOptions.saveBinaryData)->
        default_value(false),
        "Save binary data (Mode 3 assembly only).")
        ;

}



// Add configurable options to the Boost option description object.
// These are options that can be used on the command line
// and in the configuration file.
// If used in both places, the value given on the command line
// takes precedence.
void AssemblerOptions::addConfigurableOptions()
{
    using boost::program_options::value;
    using boost::program_options::bool_switch;

    configurableOptionsDescription.add_options()

        ("Reads.representation",
        value<uint64_t>(&readsOptions.representation)->
        default_value(1),
        "Read representation: 0 = raw sequence, 1 (default) = Run-Length Encoded (RLE) sequence. "
        "Experimental. Do not use.")

        ("Reads.minReadLength",
        value<int>(&readsOptions.minReadLength)->
        default_value(10000),
        "Read length cutoff. Shorter reads are discarded.")

        ("Reads.noCache",
        bool_switch(&readsOptions.noCache)->
        default_value(false),
        "If set, skip the Linux cache when loading reads. "
        "This is done by specifying the O_DIRECT flag when opening "
        "input files containing reads.")

        ("Reads.handleDuplicates",
        value<string>(&readsOptions.handleDuplicates)->
        default_value("useOneCopy"),
        "Controls handling of reads with duplicate names. "
        "Can be one of: useAllCopies, useOneCopy, useNone, forbid.")

         ("Kmers.generationMethod",
         value<int>(&kmersOptions.generationMethod)->
         default_value(0),
         "Method to generate marker k-mers: "
         "0 = random, "
         "1 = random, excluding globally overenriched k-mers,"
         "2 = random, excluding k-mers overenriched even in a single read,"
         "3 = read from file."
         "4 = random, excluding k-mers appearing in two copies close to each other even in a single read.")

         ("Kmers.k",
         value<int>(&kmersOptions.k)->
         default_value(10),
         "Length of marker k-mers (in run-length space).")

         ("Kmers.probability",
         value<double>(&kmersOptions.probability)->
         default_value(0.1, "0.1"),
         "Fraction k-mers used as a marker.")

        ("Kmers.enrichmentThreshold",
        value<double>(&kmersOptions.enrichmentThreshold)->
        default_value(100., "100."),
        "Enrichment threshold for Kmers.generationMethod 1 and 2.")

        ("Kmers.distanceThreshold",
        value<uint64_t>(&kmersOptions.distanceThreshold)->
        default_value(1000),
        "Distance threshold, in RLE bases, for Kmers.generationMethod 4")

        ("Kmers.file",
        value<string>(&kmersOptions.file),
        "The absolute path of a file containing the k-mers "
        "to be used as markers, one per line. "
        "A relative path is not accepted. "
        "Only used if Kmers.generationMethod is 3.")

        ("Kmers.globalFrequencyOverrideDirectory",
        value<string>(&kmersOptions.globalFrequencyOverrideDirectory),
        "The directory containing the hash table with marker k-mer global frequencies. "
        "Only used for Shasta development.")

        ("Assembly.crossEdgeCoverageThreshold",
        value<int>(&assemblyOptions.crossEdgeCoverageThreshold)->
        default_value(3),
        "Maximum average edge coverage for a cross edge "
        "of the assembly graph to be removed.")

        ("Assembly.markerGraphEdgeLengthThresholdForConsensus",
        value<int>(&assemblyOptions.markerGraphEdgeLengthThresholdForConsensus)->
        default_value(1000),
        "Controls assembly of long marker graph edges.")

        ("Assembly.consensusCaller",
        value<string>(&assemblyOptions.consensusCallerString)->
        default_value("Modal"),
        "Selects the consensus caller for repeat counts. "
        "See the documentation for available choices.")

        ("Assembly.storeCoverageData",
        bool_switch(&assemblyOptions.storeCoverageData)->
        default_value(false),
        "Used to request storing coverage data in binary format.")

        ("Assembly.storeCoverageDataCsvLengthThreshold",
        value<int>(&assemblyOptions.storeCoverageDataCsvLengthThreshold)->
        default_value(0),
        "Used to specify the minimum length of an assembled segment "
        "for which coverage data in csv format should be stored. "
        "If 0, no coverage data in csv format is stored.")

        ("Assembly.writeReadsByAssembledSegment",
        bool_switch(&assemblyOptions.writeReadsByAssembledSegment)->
        default_value(false),
        "Used to request writing the reads that contributed to assembling each segment.")

        ("Assembly.pruneLength",
        value<uint64_t>(&assemblyOptions.pruneLength)->
        default_value(0),
        "Prune length (in markers) for pruning of the assembly graph. "
        "Assembly graph leaves shorter than this number of markers are iteratively pruned. "
        "Set to zero to suppress pruning of the assembly graph. "
        "Assembly graph pruning takes place separately and in addition to marker graph pruning.")

        ("Assembly.detangleMethod",
        value<int>(&assemblyOptions.detangleMethod)->
        default_value(0),
        "Specify the method used to detangle the assembly graph. "
        "0 = no detangling, 1 = strict detangling, "
        "2 = less strict detangling, controlled by Assembly.detangle.* options.")

        ("Assembly.detangle.diagonalReadCountMin",
        value<uint64_t>(&assemblyOptions.detangleDiagonalReadCountMin)->
        default_value(1),
        "Minimum number of reads on detangle matrix diagonal elements "
        "required for detangling.")

        ("Assembly.detangle.offDiagonalReadCountMax",
        value<uint64_t>(&assemblyOptions.detangleOffDiagonalReadCountMax)->
        default_value(2),
        "Maximum number of reads on detangle matrix off-diagonal elements "
        "allowed for detangling.")

        ("Assembly.detangle.offDiagonalRatio",
        value<double>(&assemblyOptions.detangleOffDiagonalRatio)->
        default_value(0.3),
        "Maximum ratio of total off-diagonal elements over diagonal element "
        "allowed for detangling.")

        ("Assembly.iterative",
        bool_switch(&assemblyOptions.iterative)->
        default_value(false),
        "Used to request iterative assembly (experimental).")

        ("Assembly.iterative.iterationCount",
        value<uint64_t>(&assemblyOptions.iterativeIterationCount)->
        default_value(3),
        "Number of iterations for iterative assembly (experimental).")

        ("Assembly.iterative.pseudoPathAlignMatchScore",
        value<int64_t>(&assemblyOptions.iterativePseudoPathAlignMatchScore)->
        default_value(1),
        "Pseudopath alignment match score for iterative assembly (experimental).")

        ("Assembly.iterative.pseudoPathAlignMismatchScore",
        value<int64_t>(&assemblyOptions.iterativePseudoPathAlignMismatchScore)->
        default_value(-1),
        "Pseudopath alignment mismatch score for iterative assembly (experimental).")

        ("Assembly.iterative.pseudoPathAlignGapScore",
        value<int64_t>(&assemblyOptions.iterativePseudoPathAlignGapScore)->
        default_value(-1),
        "Pseudopath alignment gap score for iterative assembly (experimental).")

        ("Assembly.iterative.mismatchSquareFactor",
        value<double>(&assemblyOptions.iterativeMismatchSquareFactor)->
        default_value(3.),
        "Mismatch square factor for iterative assembly (experimental).")

        ("Assembly.iterative.minScore",
        value<double>(&assemblyOptions.iterativeMinScore)->
        default_value(0.),
        "Minimum pseudo-alignment score for iterative assembly (experimental).")

        ("Assembly.iterative.maxAlignmentCount",
        value<uint64_t>(&assemblyOptions.iterativeMaxAlignmentCount)->
        default_value(6),
        "Maximum number of read graph neighbors for iterative assembly (experimental).")

        ("Assembly.iterative.bridgeRemovalIterationCount",
        value<uint64_t>(&assemblyOptions.iterativeBridgeRemovalIterationCount)->
        default_value(3),
        "Number of read graph bridge removal iterations for iterative assembly (experimental).")

        ("Assembly.iterative.bridgeRemovalMaxDistance",
        value<uint64_t>(&assemblyOptions.iterativeBridgeRemovalMaxDistance)->
        default_value(2),
        "Maximum distance for read graph bridge removal for iterative assembly (experimental).")

        ("Assembly.mode3.anchorCreationMethod",
        value<string>(&assemblyOptions.mode3Options.anchorCreationMethod)->
        default_value("FromMarkerGraphEdges"),
        "Selects the method used to create anchors for mode 3 assembly. "
        "Can be: FromMarkerGraphEdges, FromMarkerKmers, FromJson.")

        ("Assembly.mode3.minAnchorCoverage",
        value<uint64_t>(&assemblyOptions.mode3Options.minAnchorCoverage)->
        default_value(0),
        "Minimum anchor coverage. "
        "If minAnchorCoverage and maxAnchorCoverage are both 0, "
        "they are set automatically to appropriate values using a simple heuristic. "
        "Only used with --Assembly.mode 3.")

        ("Assembly.mode3.maxAnchorCoverage",
        value<uint64_t>(&assemblyOptions.mode3Options.maxAnchorCoverage)->
        default_value(0),
        "Maximum anchor coverage. "
        "If minAnchorCoverage and maxAnchorCoverage are both 0, "
        "they are set automatically to appropriate values using a simple heuristic."
        "Only used with --Assembly.mode 3.")

        ("Assembly.mode3.minAnchorCoverageMultiplier",
        value<double>(&assemblyOptions.mode3Options.minAnchorCoverageMultiplier)->
        default_value(1.),
        "Multiplier applied to heuristically determined minimum anchor coverage "
        "if minAnchorCoverage and maxAnchorCoverage are both 0. "
        "Only used with --Assembly.mode 3.")

        ("Assembly.mode3.maxAnchorCoverageMultiplier",
        value<double>(&assemblyOptions.mode3Options.maxAnchorCoverageMultiplier)->
        default_value(1.),
        "Multiplier applied to heuristically determined maximum anchor coverage "
        "if minAnchorCoverage and maxAnchorCoverage are both 0. "
        "Only used with --Assembly.mode 3.")

        ("Assembly.mode3.primaryGraph.maxLoss",
        value<double>(&assemblyOptions.mode3Options.primaryGraphOptions.maxLoss)->
        default_value(0.1),
        "Use for weak edge removal in the primary graph. "
        "(Mode 3 assembly only).")

        ("Assembly.mode3.anchorGraph.crossEdgesLowCoverageThreshold",
        value<uint64_t>(&assemblyOptions.mode3Options.primaryGraphOptions.crossEdgesLowCoverageThreshold)->
        default_value(1),
        "Low coverage threshold for cross edge removal in the anchor graph. "
        "(Mode 3 assembly only).")

        ("Assembly.mode3.anchorGraph.crossEdgesHighCoverageThreshold",
        value<uint64_t>(&assemblyOptions.mode3Options.primaryGraphOptions.crossEdgesHighCoverageThreshold)->
        default_value(3),
        "High coverage threshold for cross edge removal in the anchor graph. "
        "(Mode 3 assembly only).")

        ("Assembly.mode3.assemblyGraph.detangleToleranceLow",
        value<uint64_t>(&assemblyOptions.mode3Options.assemblyGraphOptions.detangleToleranceLow)->
        default_value(0),
        "Used for detangling of the assembly graph "
        "(Mode 3 assembly only).")

        ("Assembly.mode3.assemblyGraph.detangleToleranceHigh",
        value<uint64_t>(&assemblyOptions.mode3Options.assemblyGraphOptions.detangleToleranceHigh)->
        default_value(2),
        "Used for detangling of the assembly graph "
        "(Mode 3 assembly only).")

        ("Assembly.mode3.assemblyGraph.epsilon",
        value<double>(&assemblyOptions.mode3Options.assemblyGraphOptions.epsilon)->
        default_value(0.1),
        "Epsilon value for the Bayesian model used for detangling the assembly graph "
        "(Mode 3 assembly only).")

        ("Assembly.mode3.assemblyGraph.minLogP",
        value<double>(&assemblyOptions.mode3Options.assemblyGraphOptions.minLogP)->
        default_value(20.),
        "MinLogP value (in dB) for the Bayesian model used for detangling the assembly graph "
        "(Mode 3 assembly only).")

        ("Assembly.mode3.assemblyGraph.longBubbleThreshold",
        value<uint64_t>(&assemblyOptions.mode3Options.assemblyGraphOptions.longBubbleThreshold)->
        default_value(5000),
        "Long bubble threshold "
        "(Mode 3 assembly only).")

        ("Assembly.mode3.assemblyGraph.phaseErrorThreshold",
        value<double>(&assemblyOptions.mode3Options.assemblyGraphOptions.phaseErrorThreshold)->
        default_value(0.1),
        "Phase error threshold for phasing "
        "(Mode 3 assembly only).")

        ("Assembly.mode3.assemblyGraph.bubbleErrorThreshold",
        value<double>(&assemblyOptions.mode3Options.assemblyGraphOptions.bubbleErrorThreshold)->
        default_value(0.03),
        "Bubble error threshold for bubble cleanup "
        "(Mode 3 assembly only).")

        ("Assembly.mode3.assemblyGraph.bubbleCleanupMaxOffset",
        value<uint64_t>(&assemblyOptions.mode3Options.assemblyGraphOptions.bubbleCleanupMaxOffset)->
        default_value(1000),
        "Maximum bubble offset for bubble cleanup "
        "(Mode 3 assembly only).")

        ("Assembly.mode3.assemblyGraph.chainTerminalCommonThreshold",
        value<uint64_t>(&assemblyOptions.mode3Options.assemblyGraphOptions.chainTerminalCommonThreshold)->
        default_value(3),
        "Used for bubble cleanup "
        "(Mode 3 assembly only).")

        ("Assembly.mode3.assemblyGraph.superbubbleLengthThreshold1",
        value<uint64_t>(&assemblyOptions.mode3Options.assemblyGraphOptions.superbubbleLengthThreshold1)->
        default_value(30000),
        "Length threshold used for superbubble cleanup "
        "(Mode 3 assembly only).")

        ("Assembly.mode3.assemblyGraph.superbubbleLengthThreshold2",
        value<uint64_t>(&assemblyOptions.mode3Options.assemblyGraphOptions.superbubbleLengthThreshold2)->
        default_value(10000),
        "Low length threshold used for superbubble removal "
        "(Mode 3 assembly only).")

        ("Assembly.mode3.assemblyGraph.superbubbleLengthThreshold3",
        value<uint64_t>(&assemblyOptions.mode3Options.assemblyGraphOptions.superbubbleLengthThreshold3)->
        default_value(30000),
        "High length threshold used for superbubble removal "
        "(Mode 3 assembly only).")

        ("Assembly.mode3.assemblyGraph.superbubbleLengthThreshold4",
        value<uint64_t>(&assemblyOptions.mode3Options.assemblyGraphOptions.superbubbleLengthThreshold4)->
        default_value(30000),
        "Length threshold used for superbubble detangling "
        "(Mode 3 assembly only).")

        ("Assembly.mode3.assemblyGraph.pruneLength",
        value<uint64_t>(&assemblyOptions.mode3Options.assemblyGraphOptions.pruneLength)->
        default_value(100000),
        "Length threshold used for pruning of the assembly graph."
        "(Mode 3 assembly only).")

        ("Assembly.mode3.localAssembly.estimatedOffsetRatio",
        value<double>(&assemblyOptions.mode3Options.localAssemblyOptions.estimatedOffsetRatio)->
        default_value(1.1),
        "For local assembly, the estimated offset between edgeIdA and edgeIdB gets "
        "extended by this ratio to decide how much to extend reads that only appear in edgeIdA or edgeIdB. "
        "(Mode 3 assembly only).")

        ("Assembly.mode3.localAssembly.vertexSamplingRate",
        value<double>(&assemblyOptions.mode3Options.localAssemblyOptions.vertexSamplingRate)->
        default_value(0.8),
        "Vertex sampling rate for local assembly, used to set minVertexCoverage. "
        "Only used if minVertexCoverage is 0 on input to mode3::LocalAssembly constructor. "
        "(Mode 3 assembly only).")

        ("Assembly.mode3.localAssembly.matchScore",
        value<int64_t>(&assemblyOptions.mode3Options.localAssemblyOptions.matchScore)->
        default_value(6),
        "Match score for local assembly. (Mode 3 assembly only).")

        ("Assembly.mode3.localAssembly.mismatchScore",
        value<int64_t>(&assemblyOptions.mode3Options.localAssemblyOptions.mismatchScore)->
        default_value(-1),
        "Mismatch score for local assembly. (Mode 3 assembly only).")

        ("Assembly.mode3.localAssembly.gapScore",
        value<int64_t>(&assemblyOptions.mode3Options.localAssemblyOptions.gapScore)->
        default_value(-1),
        "Gap score for local assembly. (Mode 3 assembly only).")

        ("Assembly.mode3.localAssembly.maxSkipBases",
        value<uint64_t>(&assemblyOptions.mode3Options.localAssemblyOptions.maxSkipBases)->
        default_value(500),
        "Number of bases (not markers) that can be skipped by an alignment in local assembly. "
        "(Mode 3 assembly only).")

        ("Assembly.mode3.localAssembly.maxDrift",
        value<double>(&assemblyOptions.mode3Options.localAssemblyOptions.maxDrift)->
        default_value(0.005),
        "The maximum tolerated length drift of each read. "
        "Used to compute the band for banded alignments in local assembly. "
        "(Mode 3 assembly only).")

        ("Assembly.mode3.localAssembly.minHalfBand",
        value<uint64_t>(&assemblyOptions.mode3Options.localAssemblyOptions.minHalfBand)->
        default_value(100),
        "Minimum half band, in markers, for local assembly. "
        "(Mode 3 assembly only).")

        ("Assembly.mode3.localAssembly.minScoreRatio",
        value<double>(&assemblyOptions.mode3Options.localAssemblyOptions.minScoreRatio)->
        default_value(0.7),
        "Score threshold for discarding alignments in for local assembly. "
        "(Mode 3 assembly only).")

        ("Assembly.mode3.localAssembly.maxMsaLength",
        value<uint64_t>(&assemblyOptions.mode3Options.localAssemblyOptions.maxMsaLength)->
        default_value(5000),
        "Maximum length of a multiple sequence alignment for local assembly. "
        "(Mode 3 assembly only).")

        ;
}



void ReadsOptions::write(ostream& s) const
{
    s << "[Reads]\n";
    s << "representation = " << representation << "\n";
    s << "minReadLength = " << minReadLength << "\n";
    s << "noCache = " <<
        convertBoolToPythonString(noCache) << "\n";
    s << "handleDuplicates = " << handleDuplicates << "\n";
}



void KmersOptions::write(ostream& s) const
{
    s << "[Kmers]\n";
    s << "generationMethod = " << generationMethod << "\n";
    s << "k = " << k << "\n";
    s << "probability = " << probability << "\n";
    s << "enrichmentThreshold = " << enrichmentThreshold << "\n";
    s << "distanceThreshold = " << distanceThreshold << "\n";
    s << "file = " << file << "\n";
    s << "globalFrequencyOverrideDirectory = " << globalFrequencyOverrideDirectory << "\n";
}



void AssemblyOptions::write(ostream& s) const
{
    s << "[Assembly]\n";
    s << "mode = " << mode << "\n";
    s << "crossEdgeCoverageThreshold = " << crossEdgeCoverageThreshold << "\n";
    s << "markerGraphEdgeLengthThresholdForConsensus = " <<
        markerGraphEdgeLengthThresholdForConsensus << "\n";
    s << "consensusCaller = " <<
        consensusCaller << "\n";
    s << "storeCoverageData = " <<
        convertBoolToPythonString(storeCoverageData) << "\n";
    s << "storeCoverageDataCsvLengthThreshold = " <<
        storeCoverageDataCsvLengthThreshold << "\n";
    s << "writeReadsByAssembledSegment = " <<
        convertBoolToPythonString(writeReadsByAssembledSegment) << "\n";
    s << "pruneLength = " << pruneLength << "\n";
    s << "detangleMethod = " << detangleMethod << "\n";
    s << "detangle.diagonalReadCountMin = " << detangleDiagonalReadCountMin << "\n";
    s << "detangle.offDiagonalReadCountMax = " << detangleOffDiagonalReadCountMax << "\n";
    s << "detangle.offDiagonalRatio = " << detangleOffDiagonalRatio << "\n";
    s << "iterative = " <<
        convertBoolToPythonString(iterative) << "\n";
    s << "iterative.iterationCount = " << iterativeIterationCount << "\n";
    s << "iterative.pseudoPathAlignMatchScore = " << iterativePseudoPathAlignMatchScore << "\n";
    s << "iterative.pseudoPathAlignMismatchScore = " << iterativePseudoPathAlignMismatchScore << "\n";
    s << "iterative.pseudoPathAlignGapScore = " << iterativePseudoPathAlignGapScore << "\n";
    s << "iterative.mismatchSquareFactor = " << iterativeMismatchSquareFactor << "\n";
    s << "iterative.minScore = " << iterativeMinScore << "\n";
    s << "iterative.maxAlignmentCount = " << iterativeMaxAlignmentCount << "\n";
    s << "iterative.bridgeRemovalIterationCount = " << iterativeBridgeRemovalIterationCount << "\n";
    s << "iterative.bridgeRemovalMaxDistance = " << iterativeBridgeRemovalMaxDistance << "\n";
}



void Mode3AssemblyOptions::write(ostream& s) const
{
    s << "mode3.anchorCreationMethod = " << anchorCreationMethod << "\n";
    s << "mode3.minAnchorCoverage = " << minAnchorCoverage << "\n";
    s << "mode3.maxAnchorCoverage = " << maxAnchorCoverage << "\n";
    s << "mode3.minAnchorCoverageMultiplier = " << minAnchorCoverageMultiplier << "\n";
    s << "mode3.maxAnchorCoverageMultiplier = " << maxAnchorCoverageMultiplier << "\n";
    primaryGraphOptions.write(s);
    assemblyGraphOptions.write(s);
    localAssemblyOptions.write(s);
}



void Mode3AssemblyOptions::PrimaryGraphOptions::write(ostream& s) const
{
    s << "mode3.primaryGraph.maxLoss = " << maxLoss << "\n";
    s << "mode3.anchorGraph.crossEdgesLowCoverageThreshold = " << crossEdgesLowCoverageThreshold << "\n";
    s << "mode3.anchorGraph.crossEdgesHighCoverageThreshold = " << crossEdgesHighCoverageThreshold << "\n";

}



void Mode3AssemblyOptions::AssemblyGraphOptions::write(ostream& s) const
{
    s << "mode3.assemblyGraph.detangleToleranceLow = " << detangleToleranceLow << "\n";
    s << "mode3.assemblyGraph.detangleToleranceHigh = " << detangleToleranceHigh << "\n";
    s << "mode3.assemblyGraph.epsilon = " << epsilon << "\n";
    s << "mode3.assemblyGraph.minLogP = " << minLogP << "\n";
    s << "mode3.assemblyGraph.longBubbleThreshold = " << longBubbleThreshold << "\n";
    s << "mode3.assemblyGraph.phaseErrorThreshold = " << phaseErrorThreshold << "\n";
    s << "mode3.assemblyGraph.bubbleErrorThreshold = " << bubbleErrorThreshold << "\n";
    s << "mode3.assemblyGraph.bubbleCleanupMaxOffset = " << bubbleCleanupMaxOffset << "\n";
    s << "mode3.assemblyGraph.chainTerminalCommonThreshold = " << chainTerminalCommonThreshold << "\n";
    s << "mode3.assemblyGraph.superbubbleLengthThreshold1 = " << superbubbleLengthThreshold1 << "\n";
    s << "mode3.assemblyGraph.superbubbleLengthThreshold2 = " << superbubbleLengthThreshold2 << "\n";
    s << "mode3.assemblyGraph.superbubbleLengthThreshold3 = " << superbubbleLengthThreshold3 << "\n";
    s << "mode3.assemblyGraph.superbubbleLengthThreshold4 = " << superbubbleLengthThreshold4 << "\n";
    s << "mode3.assemblyGraph.pruneLength = " << pruneLength << "\n";
}



void Mode3AssemblyOptions::LocalAssemblyOptions::write(ostream& s) const
{
    s << "mode3.localAssembly.estimatedOffsetRatio = " << estimatedOffsetRatio << "\n";
    s << "mode3.localAssembly.vertexSamplingRate = " << vertexSamplingRate << "\n";

    s << "mode3.localAssembly.matchScore = " << matchScore << "\n";
    s << "mode3.localAssembly.mismatchScore = " << mismatchScore << "\n";
    s << "mode3.localAssembly.gapScore = " << gapScore << "\n";

    s << "mode3.localAssembly.maxSkipBases = " << maxSkipBases << "\n";

    s << "mode3.localAssembly.maxDrift = " << maxDrift << "\n";

    s << "mode3.localAssembly.minHalfBand = " << minHalfBand << "\n";

    s << "mode3.localAssembly.minScoreRatio = " << minScoreRatio << "\n";

    s << "mode3.localAssembly.maxMsaLength = " << maxMsaLength << "\n";
}



void AssemblerOptions::write(ostream& s) const
{
    readsOptions.write(s);
    s << "\n";
    kmersOptions.write(s);
    s << "\n";
    assemblyOptions.write(s);
    s << endl;
}



// Function to convert a bool to True or False for better
// compatibility with Python scripts.
string shasta::convertBoolToPythonString(bool flag)
{
    return flag ? "True" : "False";
}

