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

        ("Reads.minReadLength",
        value<int>(&readsOptions.minReadLength)->
        default_value(10000),
        "Read length cutoff. Shorter reads are discarded.")

         ("Kmers.k",
         value<int>(&kmersOptions.k)->
         default_value(10),
         "Length of marker k-mers (in run-length space).")

         ("Kmers.probability",
         value<double>(&kmersOptions.probability)->
         default_value(0.1, "0.1"),
         "Fraction k-mers used as a marker.")

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

