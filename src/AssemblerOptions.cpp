#include "AssemblerOptions.hpp"
#include "ConfigurationTable.hpp"
#include "buildId.hpp"
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

    // If version was requested, write the build id and exit.
    if (variablesMap.count("version")) {
        cout << buildId() << endl;
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



    // Parse MarkerGraph.simplifyMaxLength.
    markerGraphOptions.parseSimplifyMaxLength();

    // Parse ReadOptions.desiredCoverageString into its numeric value.
    readsOptions.parseDesiredCoverageString();

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

    // Parse MarkerGraph.simplifyMaxLength.
    markerGraphOptions.parseSimplifyMaxLength();

    // Parse ReadOptions.desiredCoverageString into its numeric value.
    readsOptions.parseDesiredCoverageString();

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

        ("Reads.desiredCoverage",
        value<string>(&readsOptions.desiredCoverageString)->
        default_value("0"),
        "Reduce coverage to desired value. "
        "If not zero, specifies desired coverage (number of bases). "
        "The read length cutoff specified via --Reads.minReadLength "
        "is increased to reduce coverage to the specified value. "
        "Power of 10 multipliers can be used, for example 120Gb to "
        "request 120 Gb of coverage.")

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

        ("Reads.palindromicReads.skipFlagging",
        bool_switch(&readsOptions.palindromicReads.skipFlagging)->
        default_value(false),
        "Skip flagging palindromic reads. Oxford Nanopore reads should be flagged for better results.")

        ("Reads.palindromicReads.maxSkip",
        value<int>(&readsOptions.palindromicReads.maxSkip)->
        default_value(100),
        "Used for palindromic read detection.")

        ("Reads.palindromicReads.maxDrift",
        value<int>(&readsOptions.palindromicReads.maxDrift)->
        default_value(100),
        "Used for palindromic read detection.")

        ("Reads.palindromicReads.maxMarkerFrequency",
        value<int>(&readsOptions.palindromicReads.maxMarkerFrequency)->
        default_value(10),
        "Used for palindromic read detection.")

        ("Reads.palindromicReads.alignedFractionThreshold",
        value<double>(&readsOptions.palindromicReads.alignedFractionThreshold)->
        default_value(0.1, "0.1"),
        "Used for palindromic read detection.")

        ("Reads.palindromicReads.nearDiagonalFractionThreshold",
        value<double>(&readsOptions.palindromicReads.nearDiagonalFractionThreshold)->
        default_value(0.1, "0.1"),
        "Used for palindromic read detection.")

        ("Reads.palindromicReads.deltaThreshold",
         value<int>(&readsOptions.palindromicReads.deltaThreshold)->
         default_value(100),
         "Used for palindromic read detection.")

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

        ("MinHash.version",
        value<int>(&minHashOptions.version)->
        default_value(0),
        "Controls the version of the LowHash algorithm to use. Must be 0 (default).")

        ("MinHash.m",
        value<int>(&minHashOptions.m)->
        default_value(4),
        "The number of consecutive markers that define a MinHash/LowHash feature.")

        ("MinHash.hashFraction",
        value<double>(&minHashOptions.hashFraction)->
        default_value(0.01, "0.01"),
        "Defines how low a hash has to be to be used with the LowHash algorithm.")

        ("MinHash.minHashIterationCount",
        value<int>(&minHashOptions.minHashIterationCount)->
        default_value(10),
        "The number of MinHash/LowHash iterations, or 0 to let "
        "--MinHash.alignmentCandidatesPerRead control the number of iterations.")

        ("MinHash.alignmentCandidatesPerRead",
        value<double>(&minHashOptions.alignmentCandidatesPerRead)->
        default_value(20.),
        "If --MinHash.minHashIterationCount is 0, MinHash iteration is stopped "
        "when the average number of alignment candidates that each read is involved in "
        "reaches this value. If --MinHash.minHashIterationCount is not 0, "
        "this is not used.")

        ("MinHash.minBucketSize",
        value<int>(&minHashOptions.minBucketSize)->
        default_value(0),
        "The minimum bucket size to be used by the LowHash algorithm. "
        "If minBucketSize and maxBucketSize are both 0, they are adjusted automatically "
        "at each iteration using simple heuristics.")

        ("MinHash.maxBucketSize",
        value<int>(&minHashOptions.maxBucketSize)->
        default_value(10),
        "The maximum bucket size to be used by the LowHash algorithm. "
        "If minBucketSize and maxBucketSize are both 0, they are adjusted automatically "
        "at each iteration using simple heuristics.")

        ("MinHash.minFrequency",
        value<int>(&minHashOptions.minFrequency)->
        default_value(2),
        "The minimum number of times a pair of reads must be found by the MinHash/LowHash algorithm "
        "in order to be considered a candidate alignment.")

        ("MinHash.allPairs",
        bool_switch(&minHashOptions.allPairs)->
        default_value(false),
        "Skip the MinHash algorithm and mark all pairs of reads as alignment"
        "candidates with both orientation. This should only be used for experimentation "
        "on very small runs because it is very time consuming.")

        ("Align.alignMethod",
        value<int>(&alignOptions.alignMethod)->
        default_value(3),
        "The alignment method to be used to create the read graph & the marker graph. "
        "0 = old Shasta method, 1 = SeqAn (slow), 3 = banded SeqAn, 4, 5, 6 = experimental.")

        ("Align.maxSkip",
        value<int>(&alignOptions.maxSkip)->
        default_value(30),
        "The maximum number of markers that an alignment is allowed to skip.")

        ("Align.maxDrift",
        value<int>(&alignOptions.maxDrift)->
        default_value(30),
        "The maximum amount of marker drift that an alignment is allowed to tolerate "
        "between successive markers.")

        ("Align.maxTrim",
        value<int>(&alignOptions.maxTrim)->
        default_value(30),
        "The maximum number of unaligned markers tolerated at the beginning and end of an alignment.")

        ("Align.maxMarkerFrequency",
        value<int>(&alignOptions.maxMarkerFrequency)->
        default_value(10),
        "Marker frequency threshold. Markers more frequent than this value in either of "
        "two oriented reads being aligned are discarded and not used to compute "
        "the alignment.")

        ("Align.minAlignedMarkerCount",
        value<int>(&alignOptions.minAlignedMarkerCount)->
        default_value(100),
        "The minimum number of aligned markers for an alignment to be used.")

        ("Align.minAlignedFraction",
        value<double>(&alignOptions.minAlignedFraction)->
        default_value(0.),
        "The minimum fraction of aligned markers for an alignment to be used.")

        ("Align.matchScore",
        value<int>(&alignOptions.matchScore)->
        default_value(6),
        "Match score for marker alignments (only used for alignment methods 1 and 3).")

        ("Align.mismatchScore",
        value<int>(&alignOptions.mismatchScore)->
        default_value(-1),
        "Mismatch score for marker alignments (only used for alignment methods 1 and 3).")

        ("Align.gapScore",
        value<int>(&alignOptions.gapScore)->
        default_value(-1),
        "Gap score for marker alignments (only used for alignment methods 1 and 3).")

        ("Align.downsamplingFactor",
        value<double>(&alignOptions.downsamplingFactor)->
        default_value(0.1),
        "Downsampling factor (only used for alignment method 3).")

        ("Align.bandExtend",
        value<int>(&alignOptions.bandExtend)->
        default_value(10),
        "Amount to extend the downsampled band "
        "(only used for alignment method 3).")

        ("Align.maxBand",
        value<int>(&alignOptions.maxBand)->
        default_value(1000),
        "Maximum alignment band "
        "(only used for alignment method 3).")

        ("Align.sameChannelReadAlignment.suppressDeltaThreshold",
        value<int>(&alignOptions.sameChannelReadAlignmentSuppressDeltaThreshold)->
        default_value(0),
        "If not zero, alignments between reads from the same nanopore channel "
        "and close in time are suppressed. The \"read\" meta data fields "
        "from the FASTA or FASTQ header are checked. If their difference, in "
        "absolute value, is less than the value of this option, the alignment "
        "is suppressed. This can help avoid assembly artifact. "
        "This check is only done if the two reads have identical meta data fields "
        "\"runid\", \"sampleid\", and \"ch\". "
        "If any of these meta data fields are missing, this check is suppressed and this "
        "option has no effect.")

        ("Align.suppressContainments",
        bool_switch(&alignOptions.suppressContainments)->
        default_value(false),
        "Suppress containment alignments, that is alignments in which "
        "one read is entirely contained in another read, "
        "except possibly for up to maxTrim markers at the beginning and end.")

        ("Align.align4.deltaX",
        value<uint64_t>(&alignOptions.align4DeltaX)->
        default_value(200),
        "Only used for alignment method 4 (experimental).")

        ("Align.align4.deltaY",
        value<uint64_t>(&alignOptions.align4DeltaY)->
        default_value(10),
        "Only used for alignment method 4 (experimental).")

        ("Align.align4.minEntryCountPerCell",
        value<uint64_t>(&alignOptions.align4MinEntryCountPerCell)->
        default_value(10),
        "Only used for alignment method 4 (experimental).")

        ("Align.align4.maxDistanceFromBoundary",
        value<uint64_t>(&alignOptions.align4MaxDistanceFromBoundary)->
        default_value(100),
        "Only used for alignment method 4 (experimental).")

        ("Align.align5.driftRateTolerance",
        value<double>(&alignOptions.align5DriftRateTolerance)->
        default_value(0.02),
        "Maximum allowed drift rate for alignment method 5.")

        ("Align.align5.minBandExtend",
        value<uint64_t>(&alignOptions.align5MinBandExtend)->
        default_value(10),
        "Minimum band extension for alignment method 5.")

        ("Align.align6.maxLocalFrequency",
        value<uint64_t>(&alignOptions.align6Options.maxLocalFrequency)->
        default_value(100),
        "Only used for alignment method 6.")

        ("Align.align6.minGlobalFrequency",
        value<uint64_t>(&alignOptions.align6Options.minGlobalFrequency)->
        default_value(0),
        "Only used for alignment method 6.")

        ("Align.align6.maxGlobalFrequency",
        value<uint64_t>(&alignOptions.align6Options.maxGlobalFrequency)->
        default_value(0),
        "Only used for alignment method 6.")

        ("Align.align6.maxGlobalFrequencyMultiplier",
        value<double>(&alignOptions.align6Options.maxGlobalFrequencyMultiplier)->
        default_value(2.),
        "Only used for alignment method 6.")

        ("Align.align6.minLowFrequencyCount",
        value<uint64_t>(&alignOptions.align6Options.minLowFrequencyCount)->
        default_value(4),
        "Only used for alignment method 6.")

        ("Align.align6.driftRateTolerance",
        value<double>(&alignOptions.align6Options.driftRateTolerance)->
        default_value(0.05),
        "Only used for alignment method 6.")

        ("Align.align6.maxInBandCount",
        value<uint64_t>(&alignOptions.align6Options.maxInBandCount)->
        default_value(1000),
        "Only used for alignment method 6.")

        ("Align.align6.maxInBandRatio",
        value<double>(&alignOptions.align6Options.maxInBandRatio)->
        default_value(100.),
        "Only used for alignment method 6.")

        ("ReadGraph.creationMethod",
        value<int>(&readGraphOptions.creationMethod)->
        default_value(0),
        "The method used to create the read graph (0 default, 1 or 2 experimental).")

        ("ReadGraph.maxAlignmentCount",
        value<int>(&readGraphOptions.maxAlignmentCount)->
        default_value(6),
        "The maximum number of alignments to be kept for each read.")

        ("ReadGraph.preferAlignedFraction",
        bool_switch(&readGraphOptions.preferAlignedFraction)->
        default_value(false),
        "Prefer alignments with a higher aligned fraction, rather than higher number of aligned markers.")

        ("ReadGraph.maxChimericReadDistance",
        value<int>(&readGraphOptions.maxChimericReadDistance)->
        default_value(2),
        "Used for chimeric read detection. Set to 0 to turn off chimera detection.")

        ("ReadGraph.strandSeparationMethod",
        value<uint64_t>(&readGraphOptions.strandSeparationMethod)->
        default_value(1),
        "Strand separation method: "
        "0 = no strand separation, "
        "1 = limited strand separation, "
        "2 = strict strand separation.")

        ("ReadGraph.crossStrandMaxDistance",
        value<int>(&readGraphOptions.crossStrandMaxDistance)->
        default_value(6),
        "Maximum distance (edges) for strand separation method 1.")

        ("ReadGraph.removeConflicts",
        bool_switch(&readGraphOptions.removeConflicts)->
        default_value(false),
        "Remove conflicts from the read graph. Experimental - do not use.")

        ("ReadGraph.markerCountPercentile",
        value<double>(&readGraphOptions.markerCountPercentile)->
        default_value(0.015, "0.015"),
        "Percentile for --ReadGraph.markerCount "
        "(only used when --ReadGraph.creationMethod is 2).")

        ("ReadGraph.alignedFractionPercentile",
        value<double>(&readGraphOptions.alignedFractionPercentile)->
        default_value(0.12, "0.12"),
        "Percentile for adaptive selection of --ReadGraph.alignedFraction "
        "(only used when --ReadGraph.creationMethod is 2).")

        ("ReadGraph.maxSkipPercentile",
        value<double>(&readGraphOptions.maxSkipPercentile)->
        default_value(0.12, "0.12"),
        "Percentile for adaptive selection of --ReadGraph.maxSkip "
        "(only used when --ReadGraph.creationMethod is 2).")

        ("ReadGraph.maxDriftPercentile",
        value<double>(&readGraphOptions.maxDriftPercentile)->
        default_value(0.12, "0.12"),
        "Percentile for adaptive selection of --ReadGraph.maxDrift "
        "(only used when --ReadGraph.creationMethod is 2).")

        ("ReadGraph.maxTrimPercentile",
        value<double>(&readGraphOptions.maxTrimPercentile)->
        default_value(0.015, "0.015"),
        "Percentile for adaptive selection of --ReadGraph.maxTrim "
        "(only used when --ReadGraph.creationMethod is 2).")

        ("ReadGraph.flagInconsistentAlignments",
        bool_switch(&readGraphOptions.flagInconsistentAlignments)->
        default_value(false),
        "Flag inconsistent alignments. Experimental.")

        ("ReadGraph.flagInconsistentAlignments.triangleErrorThreshold",
        value<uint64_t>(&readGraphOptions.flagInconsistentAlignmentsTriangleErrorThreshold)->
        default_value(200),
        "Triangle error threshold, in markers, for flagging inconsistent alignments. "
        "Only used if --ReadGraph.flagInconsistentAlignments is set. Experimental.")

        ("ReadGraph.flagInconsistentAlignments.leastSquareErrorThreshold",
        value<uint64_t>(&readGraphOptions.flagInconsistentAlignmentsLeastSquareErrorThreshold)->
        default_value(200),
        "Least square error threshold, in markers, for flagging inconsistent alignments. "
        "Only used if --ReadGraph.flagInconsistentAlignments is set. Experimental.")

        ("ReadGraph.flagInconsistentAlignments.leastSquareMaxDistance",
        value<uint64_t>(&readGraphOptions.flagInconsistentAlignmentsLeastSquareMaxDistance)->
        default_value(1),
        "Least square max distance for flagging inconsistent alignments. "
        "Only used if --ReadGraph.flagInconsistentAlignments is set. Experimental.")

        //* New readGraph4withStrandSeparation options

        ("ReadGraph.epsilon",
        value<double>(&readGraphOptions.epsilon)->
        default_value(1e-4, "1e-4"),
        "Mismatch rate used for the Bayesian ranking of alignments. "
        "(only used when --ReadGraph.creationMethod is 4).")

        ("ReadGraph.delta",
        value<double>(&readGraphOptions.delta)->
        default_value(5e-4, "5e-4"),
        "Divergence rate δ (fraction of discordant bases) used for the Bayesian ranking of alignments. "
        "(only used when --ReadGraph.creationMethod is 4).")

        ("ReadGraph.WThreshold",
        value<double>(&readGraphOptions.WThreshold)->
        default_value(1e-8, "1e-8"),
        "Logarithm of probability ratio used for the Bayesian ranking of alignments. "
        "(only used when --ReadGraph.creationMethod is 4).")

        ("ReadGraph.WThresholdForBreaks",
        value<double>(&readGraphOptions.WThresholdForBreaks)->
        default_value(1e+15, "1e+15"),
        "Logarithm of probability ratio used for the Bayesian ranking of alignments to detect breaks in the read graph. "
        "(only used when --ReadGraph.creationMethod is 4).")

        ("MarkerGraph.minCoverage",
        value<int>(&markerGraphOptions.minCoverage)->
        default_value(10),
        "Minimum coverage (number of supporting oriented reads) "
        "for a marker graph vertex to be created."
        "Specifying 0 causes a suitable value of this parameter "
        "to be selected automatically.")

        ("MarkerGraph.maxCoverage",
        value<int>(&markerGraphOptions.maxCoverage)->
        default_value(100),
        "Maximum coverage (number of supporting oriented reads) "
        "for a marker graph vertex.")

        ("MarkerGraph.minCoveragePerStrand",
        value<int>(&markerGraphOptions.minCoveragePerStrand)->
        default_value(0),
        "Minimum coverage (number of supporting oriented reads) "
        "for each strand for a marker graph vertex.")

        ("MarkerGraph.minEdgeCoverage",
        value<uint64_t>(&markerGraphOptions.minEdgeCoverage)->
        default_value(6),
        "Minimum edge coverage (number of supporting oriented reads) "
        "for a marker graph edge to be created."
        "Only used with --Assembly.mode 2.")

        ("MarkerGraph.minEdgeCoveragePerStrand",
        value<uint64_t>(&markerGraphOptions.minEdgeCoveragePerStrand)->
        default_value(2),
        "Minimum edge coverage (number of supporting oriented reads) "
        "on each strand "
        "for a marker graph edge to be created."
        "Only used with --Assembly.mode 2.")

        ("MarkerGraph.allowDuplicateMarkers",
        bool_switch(&markerGraphOptions.allowDuplicateMarkers)->
        default_value(false),
        "Specifies whether to allow more than one marker on the "
        "same oriented read on a single marker graph vertex. Experimental.")

        ("MarkerGraph.cleanupDuplicateMarkers",
        bool_switch(&markerGraphOptions.cleanupDuplicateMarkers)->
        default_value(false),
        "Specifies whether to clean up marker graph vertices with more than one marker on the "
        "same oriented read. Experimental.")

        ("MarkerGraph.duplicateMarkersPattern1Threshold",
        value<double>(&markerGraphOptions.duplicateMarkersPattern1Threshold)->
        default_value(0.5),
        "Used when cleaning up marker graph vertices with more than one marker on the "
        "same oriented read. Experimental.")

        ("MarkerGraph.lowCoverageThreshold",
        value<int>(&markerGraphOptions.lowCoverageThreshold)->
        default_value(0),
        "Used during approximate transitive reduction. Marker graph edges with coverage "
        "lower than this value are always marked as removed regardless of reachability.")

        ("MarkerGraph.highCoverageThreshold",
        value<int>(&markerGraphOptions.highCoverageThreshold)->
        default_value(256),
        "Used during approximate transitive reduction. Marker graph edges with coverage "
        "higher than this value are never marked as removed regardless of reachability.")

        ("MarkerGraph.maxDistance",
        value<int>(&markerGraphOptions.maxDistance)->
        default_value(30),
        "Used during approximate transitive reduction.")

        ("MarkerGraph.edgeMarkerSkipThreshold",
        value<int>(&markerGraphOptions.edgeMarkerSkipThreshold)->
        default_value(100),
        "Used during approximate transitive reduction.")

        ("MarkerGraph.pruneIterationCount",
        value<int>(&markerGraphOptions.pruneIterationCount)->
        default_value(6),
        "Number of prune iterations.")

        ("MarkerGraph.simplifyMaxLength",
        value<string>(&markerGraphOptions.simplifyMaxLength)->
        default_value("10,100,1000"),
        "Maximum lengths (in markers) used at each iteration of simplifyMarkerGraph.")

        ("MarkerGraph.crossEdgeCoverageThreshold",
        value<double>(&markerGraphOptions.crossEdgeCoverageThreshold)->
        default_value(0.),
        "Experimental. Cross edge coverage threshold. If this is not zero, assembly graph cross-edges "
        "with average edge coverage less than this value are removed, together with the "
        "corresponding marker graph edges. A cross edge is defined as an edge v0->v1 "
        "with out-degree(v0)>1, in-degree(v1)>1.")

        ("MarkerGraph.peakFinder.minAreaFraction",
        value<double>(&markerGraphOptions.peakFinderMinAreaFraction)->
        default_value(0.08),
        "Used in the automatic selection of --MarkerGraph.minCoverage when "
        "--MarkerGraph.minCoverage is set to 0.")

        ("MarkerGraph.peakFinder.areaStartIndex",
        value<uint64_t>(&markerGraphOptions.peakFinderAreaStartIndex)->
        default_value(2),
        "Used in the automatic selection of --MarkerGraph.minCoverage when "
        "--MarkerGraph.minCoverage is set to 0.")

        ("MarkerGraph.secondaryEdges.maxSkip",
        value<uint64_t>(&markerGraphOptions.secondaryEdgesMaxSkip)->
        default_value(1000000),
        "Maximum number of markers skipped by a secondary edge (mode 2 assembly only).")

        ("MarkerGraph.secondaryEdges.split.errorRateThreshold",
        value<double>(&markerGraphOptions.secondaryEdgesSplitErrorRateThreshold)->
        default_value(0.25),
        "Error rate threshold used for splitting secondary edges (mode 2 assembly only).")

        ("MarkerGraph.secondaryEdges.split.minCoverage",
        value<uint64_t>(&markerGraphOptions.secondaryEdgesSplitMinCoverage)->
        default_value(4),
        "Minimum coverage for secondary edges generated during splitting (mode 2 assembly only).")

        ("Assembly.mode",
        value<uint64_t>(&assemblyOptions.mode)->
        default_value(0),
        "Assembly mode (0=default, 1=experimental).")

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

        ("Assembly.mode2.strongBranchThreshold",
        value<uint64_t>(&assemblyOptions.mode2Options.strongBranchThreshold)->
        default_value(2),
        "Minimum number of supporting reads required for a strong branch. Only used in Mode 2 assembly.")

        ("Assembly.mode2.epsilon",
        value<double>(&assemblyOptions.mode2Options.epsilon)->
        default_value(0.1),
        "Epsilon for the Bayesian model used for phasing and for bubble removal."
        "This is the probability that a read appears on the wrong branch of a diploid bubble. "
        "Only used in Mode 2 assembly.")

        ("Assembly.mode2.bubbleRemoval.minConcordantReadCount",
        value<uint64_t>(&assemblyOptions.mode2Options.minConcordantReadCountForBubbleRemoval)->
        default_value(3),
        "Minimum number of concordant reads for bubble removal. "
        "Only used in Mode 2 assembly.")

        ("Assembly.mode2.bubbleRemoval.maxDiscordantReadCount",
        value<uint64_t>(&assemblyOptions.mode2Options.maxDiscordantReadCountForBubbleRemoval)->
        default_value(6),
        "Maximum number of discordant reads for bubble removal. "
        "Only used in Mode 2 assembly.")

        ("Assembly.mode2.bubbleRemoval.minLogP",
        value<double>(&assemblyOptions.mode2Options.minLogPForBubbleRemoval)->
        default_value(30.),
        "Minimul log(P) (in decibels) for bubble removal. "
        "Only used in Mode 2 assembly.")

        ("Assembly.mode2.bubbleRemoval.componentSizeThreshold",
        value<uint64_t>(&assemblyOptions.mode2Options.componentSizeThresholdForBubbleRemoval)->
        default_value(10),
        "Component size threshold for bubble removal. "
        "Only used in Mode 2 assembly.")

        ("Assembly.mode2.phasing.minConcordantReadCount",
        value<uint64_t>(&assemblyOptions.mode2Options.minConcordantReadCountForPhasing)->
        default_value(2),
        "Minimum number of concordant reads for phasing. "
        "Only used in Mode 2 assembly.")

        ("Assembly.mode2.phasing.maxDiscordantReadCount",
        value<uint64_t>(&assemblyOptions.mode2Options.maxDiscordantReadCountForPhasing)->
        default_value(1),
        "Maximum number of discordant reads for phasing. "
        "Only used in Mode 2 assembly.")

        ("Assembly.mode2.phasing.minLogP",
        value<double>(&assemblyOptions.mode2Options.minLogPForPhasing)->
        default_value(10.),
        "Minimul log(P) (in decibels) for phasing. "
        "Only used in Mode 2 assembly.")

        ("Assembly.mode2.superbubble.maxSize",
        value<uint64_t>(&assemblyOptions.mode2Options.maxSuperbubbleSize)->
        default_value(50),
        "Maximum size (number of edges) of a superbubble to be processed. "
        "Only used in Mode 2 assembly.")

        ("Assembly.mode2.superbubble.maxChunkSize",
        value<uint64_t>(&assemblyOptions.mode2Options.maxSuperbubbleChunkSize)->
        default_value(20),
        "Maximum size (numbef of edges) of a superbubble chunk to be processed. "
        "Only used in Mode 2 assembly.")

        ("Assembly.mode2.superbubble.maxChunkPathCount",
        value<uint64_t>(&assemblyOptions.mode2Options.maxSuperbubbleChunkPathCount)->
        default_value(20),
        "Maximum number of paths to be processed in a superbubble chunk. "
        "Only used in Mode 2 assembly.")

        ("Assembly.mode2.superbubble.edgeLengthThreshold",
        value<uint64_t>(&assemblyOptions.mode2Options.superbubbleEdgeLengthThreshold)->
        default_value(6),
        "Edge length threshold (in markers) for superbubble removal. "
        "Only used in Mode 2 assembly.")

        ("Assembly.mode2.suppressGfaOutput",
        bool_switch(&assemblyOptions.mode2Options.suppressGfaOutput)->
        default_value(false),
        "Suppress all GFA output (Mode 2 assembly only).")

        ("Assembly.mode2.suppressFastaOutput",
        bool_switch(&assemblyOptions.mode2Options.suppressFastaOutput)->
        default_value(false),
        "Suppress all FASTA output (Mode 2 assembly only).")

        ("Assembly.mode2.suppressDetailedOutput",
        bool_switch(&assemblyOptions.mode2Options.suppressDetailedOutput)->
        default_value(false),
        "Suppress output of detailed representation of the assembly (Mode 2 assembly only).")

        ("Assembly.mode2.suppressPhasedOutput",
        bool_switch(&assemblyOptions.mode2Options.suppressPhasedOutput)->
        default_value(false),
        "Suppress output of phased representation of the assembly (Mode 2 assembly only).")

        ("Assembly.mode2.suppressHaploidOutput",
        bool_switch(&assemblyOptions.mode2Options.suppressHaploidOutput)->
        default_value(false),
        "Suppress output of haploid representation of the assembly (Mode 2 assembly only).")

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



void PalindromicReadOptions::write(ostream& s) const
{
    s << "palindromicReads.skipFlagging = " << convertBoolToPythonString(skipFlagging) << "\n";
    s << "palindromicReads.maxSkip = " << maxSkip << "\n";
    s << "palindromicReads.maxDrift = " << maxDrift << "\n";
    s << "palindromicReads.maxMarkerFrequency = " << maxMarkerFrequency << "\n";
    s << "palindromicReads.alignedFractionThreshold = " << alignedFractionThreshold << "\n";
    s << "palindromicReads.nearDiagonalFractionThreshold = " << nearDiagonalFractionThreshold << "\n";
    s << "palindromicReads.deltaThreshold = " << deltaThreshold << "\n";
}



void ReadsOptions::write(ostream& s) const
{
    s << "[Reads]\n";
    s << "representation = " << representation << "\n";
    s << "minReadLength = " << minReadLength << "\n";
    s << "desiredCoverage = " << desiredCoverageString << "\n";
    s << "noCache = " <<
        convertBoolToPythonString(noCache) << "\n";
    s << "handleDuplicates = " << handleDuplicates << "\n";
    palindromicReads.write(s);
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



void MinHashOptions::write(ostream& s) const
{
    s << "[MinHash]\n";
    s << "version = " << version << "\n";
    s << "m = " << m << "\n";
    s << "hashFraction = " << hashFraction << "\n";
    s << "minHashIterationCount = " << minHashIterationCount << "\n";
    s << "alignmentCandidatesPerRead = " << alignmentCandidatesPerRead << "\n";
    s << "minBucketSize = " << minBucketSize << "\n";
    s << "maxBucketSize = " << maxBucketSize << "\n";
    s << "minFrequency = " << minFrequency << "\n";
    s << "allPairs = " <<
        convertBoolToPythonString(allPairs) << "\n";
}



void AlignOptions::write(ostream& s) const
{
    s << "[Align]\n";
    s << "alignMethod = " << alignMethod << "\n";
    s << "maxSkip = " << maxSkip << "\n";
    s << "maxDrift = " << maxDrift << "\n";
    s << "maxTrim = " << maxTrim << "\n";
    s << "maxMarkerFrequency = " << maxMarkerFrequency << "\n";
    s << "minAlignedMarkerCount = " << minAlignedMarkerCount << "\n";
    s << "minAlignedFraction = " << minAlignedFraction << "\n";
    s << "matchScore = " << matchScore << "\n";
    s << "mismatchScore = " << mismatchScore << "\n";
    s << "gapScore = " << gapScore << "\n";
    s << "downsamplingFactor = " << downsamplingFactor << "\n";
    s << "bandExtend = " << bandExtend << "\n";
    s << "maxBand = " << maxBand << "\n";
    s << "sameChannelReadAlignment.suppressDeltaThreshold = " <<
        sameChannelReadAlignmentSuppressDeltaThreshold << "\n";
    s << "suppressContainments = " <<
        convertBoolToPythonString(suppressContainments) << "\n";
    s << "align4.deltaX = " << align4DeltaX << "\n";
    s << "align4.deltaY = " << align4DeltaY << "\n";
    s << "align4.minEntryCountPerCell = " << align4MinEntryCountPerCell << "\n";
    s << "align4.maxDistanceFromBoundary = " << align4MaxDistanceFromBoundary << "\n";
    s << "align5.driftRateTolerance = " << align5DriftRateTolerance << "\n";
    s << "align5.minBandExtend = " << align5MinBandExtend << "\n";

    align6Options.write(s);

}



void Align6Options::write(ostream& s) const
{
    s << "align6.maxLocalFrequency = " <<  maxLocalFrequency << "\n";
    s << "align6.minGlobalFrequency = " << minGlobalFrequency << "\n";
    s << "align6.maxGlobalFrequency = " << maxGlobalFrequency << "\n";
    s << "align6.maxGlobalFrequencyMultiplier = " << maxGlobalFrequencyMultiplier << "\n";
    s << "align6.minLowFrequencyCount = " << minLowFrequencyCount << "\n";
    s << "align6.driftRateTolerance = " << driftRateTolerance << "\n";
    s << "align6.maxInBandCount = " << maxInBandCount << "\n";
    s << "align6.maxInBandRatio = " << maxInBandRatio << "\n";
}



void ReadGraphOptions::write(ostream& s) const
{
    s << "[ReadGraph]\n";
    s << "creationMethod = " << creationMethod << "\n";
    s << "maxAlignmentCount = " << maxAlignmentCount << "\n";
    s << "preferAlignedFraction = " << convertBoolToPythonString(preferAlignedFraction) << "\n";
    s << "maxChimericReadDistance = " << maxChimericReadDistance << "\n";
    s << "strandSeparationMethod = " << strandSeparationMethod << "\n";
    s << "crossStrandMaxDistance = " << crossStrandMaxDistance << "\n";
    s << "removeConflicts = " <<
        convertBoolToPythonString(removeConflicts) << "\n";
    s << "markerCountPercentile = " << markerCountPercentile << "\n";
    s << "alignedFractionPercentile = " << alignedFractionPercentile << "\n";
    s << "maxSkipPercentile = " << maxSkipPercentile << "\n";
    s << "maxDriftPercentile = " << maxDriftPercentile << "\n";
    s << "maxTrimPercentile = " << maxTrimPercentile << "\n";

    s << "flagInconsistentAlignments = " << convertBoolToPythonString(flagInconsistentAlignments) << "\n";
    s << "flagInconsistentAlignments.triangleErrorThreshold = " <<
        flagInconsistentAlignmentsTriangleErrorThreshold << "\n";
    s << "flagInconsistentAlignments.leastSquareErrorThreshold = " <<
        flagInconsistentAlignmentsLeastSquareErrorThreshold << "\n";
    s << "flagInconsistentAlignments.leastSquareMaxDistance = " <<
        flagInconsistentAlignmentsLeastSquareMaxDistance << "\n";
    // New readGraph4withStrandSeparation options
    s << "epsilon = " << epsilon << "\n";
    s << "delta = " << delta << "\n";
    s << "WThreshold = " << WThreshold << "\n";
    s << "WThresholdForBreaks = " << WThresholdForBreaks << "\n";
}



void MarkerGraphOptions::write(ostream& s) const
{
    s << "[MarkerGraph]\n";
    s << "minCoverage = " << minCoverage << "\n";
    s << "maxCoverage = " << maxCoverage << "\n";
    s << "minCoveragePerStrand = " << minCoveragePerStrand << "\n";
    s << "minEdgeCoverage = " << minEdgeCoverage << "\n";
    s << "minEdgeCoveragePerStrand = " << minEdgeCoveragePerStrand << "\n";
    s << "allowDuplicateMarkers = " <<
        convertBoolToPythonString(allowDuplicateMarkers) << "\n";
    s << "cleanupDuplicateMarkers = " <<
        convertBoolToPythonString(cleanupDuplicateMarkers) << "\n";
    s << "duplicateMarkersPattern1Threshold = " << duplicateMarkersPattern1Threshold << "\n";
    s << "lowCoverageThreshold = " << lowCoverageThreshold << "\n";
    s << "highCoverageThreshold = " << highCoverageThreshold << "\n";
    s << "maxDistance = " << maxDistance << "\n";
    s << "edgeMarkerSkipThreshold = " << edgeMarkerSkipThreshold << "\n";
    s << "pruneIterationCount = " << pruneIterationCount << "\n";
    s << "simplifyMaxLength = " << simplifyMaxLength << "\n";
    s << "crossEdgeCoverageThreshold = " << crossEdgeCoverageThreshold << "\n";
    s << "peakFinder.minAreaFraction = " << peakFinderMinAreaFraction << "\n";
    s << "peakFinder.areaStartIndex = " << peakFinderAreaStartIndex << "\n";

    s << "secondaryEdges.maxSkip = " << secondaryEdgesMaxSkip << "\n";
    s << "secondaryEdges.split.errorRateThreshold = " << secondaryEdgesSplitErrorRateThreshold << "\n";
    s << "secondaryEdges.split.minCoverage = " << secondaryEdgesSplitMinCoverage << "\n";
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

    mode2Options.write(s);
    mode3Options.write(s);
}



void Mode2AssemblyOptions::write(ostream& s) const
{
    s << "mode2.strongBranchThreshold = " << strongBranchThreshold << "\n";
    s << "mode2.epsilon = " << epsilon << "\n";
    s << "mode2.bubbleRemoval.minConcordantReadCount = " << minConcordantReadCountForBubbleRemoval << "\n";
    s << "mode2.bubbleRemoval.maxDiscordantReadCount = " << maxDiscordantReadCountForBubbleRemoval << "\n";
    s << "mode2.bubbleRemoval.minLogP = " << minLogPForBubbleRemoval << "\n";
    s << "mode2.bubbleRemoval.componentSizeThreshold = " << componentSizeThresholdForBubbleRemoval << "\n";
    s << "mode2.phasing.minConcordantReadCount = " << minConcordantReadCountForPhasing << "\n";
    s << "mode2.phasing.maxDiscordantReadCount = " << maxDiscordantReadCountForPhasing << "\n";
    s << "mode2.phasing.minLogP = " << minLogPForPhasing << "\n";
    s << "mode2.superbubble.maxSize = " << maxSuperbubbleSize << "\n";
    s << "mode2.superbubble.maxChunkSize = " << maxSuperbubbleChunkSize << "\n";
    s << "mode2.superbubble.maxChunkPathCount = " << maxSuperbubbleChunkPathCount << "\n";
    s << "mode2.superbubble.edgeLengthThreshold = " << superbubbleEdgeLengthThreshold << "\n";

    s << "mode2.suppressGfaOutput = " << convertBoolToPythonString(suppressGfaOutput) << "\n";
    s << "mode2.suppressFastaOutput = " << convertBoolToPythonString(suppressFastaOutput) << "\n";
    s << "mode2.suppressDetailedOutput = " << convertBoolToPythonString(suppressDetailedOutput) << "\n";
    s << "mode2.suppressPhasedOutput = " << convertBoolToPythonString(suppressPhasedOutput) << "\n";
    s << "mode2.suppressHaploidOutput = " << convertBoolToPythonString(suppressHaploidOutput) << "\n";
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
    minHashOptions.write(s);
    s << "\n";
    alignOptions.write(s);
    s << "\n";
    readGraphOptions.write(s);
    s << "\n";
    markerGraphOptions.write(s);
    s << "\n";
    assemblyOptions.write(s);
    s << endl;
}



void MarkerGraphOptions::parseSimplifyMaxLength()
{
    simplifyMaxLengthVector.clear();

    boost::tokenizer< boost::char_separator<char> > tokenizer(
        simplifyMaxLength, boost::char_separator<char>(","));
    for(const string& token: tokenizer) {
        try {
            size_t numberEndsHere;
            const size_t value = std::stoi(token, &numberEndsHere);
            if(numberEndsHere != token.size()) {
                throw runtime_error("Error parsing MarkerGraph.simplifyMaxLength " +
                    simplifyMaxLength);
            }
            simplifyMaxLengthVector.push_back(value);
        } catch(const std::invalid_argument& e) {
            throw runtime_error("Error parsing MarkerGraph,simplifyMaxLength " +
                simplifyMaxLength);
        }
    }

}



void ReadsOptions::parseDesiredCoverageString() {
    size_t pos = 0;
    desiredCoverage = std::stoull(desiredCoverageString, &pos);
    
    if (pos == desiredCoverageString.size()) {
        // Raw number of bases specified.
        return;
    }
    
    string suffix = desiredCoverageString.substr(pos);
    if (suffix == "Gbp" or suffix == "Gb" or suffix == "G") {
        desiredCoverage *= 1000 * 1000 * 1000;
    } else if (suffix == "Mbp" or suffix == "Mb" or suffix == "M") {
        desiredCoverage *= 1000 * 1000;
    } else if (suffix == "Kbp" or suffix == "Kb" or suffix == "K") {
        desiredCoverage *= 1000;
    } else {
        throw runtime_error("Unsupported units used for specifying desiredCoverage.");
    }
}



// Function to convert a bool to True or False for better
// compatibility with Python scripts.
string shasta::convertBoolToPythonString(bool flag)
{
    return flag ? "True" : "False";
}

