#ifndef SHASTA_ASSEMBLER_OPTIONS_HPP
#define SHASTA_ASSEMBLER_OPTIONS_HPP


/*******************************************************************************

CLASSES DESCRIBING ASSEMBLER OPTIONS

There are two types of options:
- Options that can be used both on the command line and in a configuration file
  ("configurable options").
- Options that can only be used in a configuration file
  ("non-configurable options").

Configuration files are divided in sections formatted like this:

[SectionName]
optionName = optionValue

The command line syntax corresponding to the above is
--SectionName.optionName optionValue

Each SectionName corresponds to a class defined below. For example,
section [Align] corresonds to class AlignOptions.

If the option is a Boolean switch, use True or False as the optionValue.



ADDING A NEW CONFIGURABLE OPTION

1. Add the option to the class corresponding to the desired section.
2. Modify the write function to that class to also write the newly added option.
3. Modify AssemblerOptions::addConfigurableOptions to reflect the new option,
   making sure to include a default value and at least a minimal help message.
4. Document the option in shasta/docs/CommandLineOptions.html.
5. If the option requires validation add it at the appropriate place in
   shasta/srcMain/main.cpp.

For options not ready for end users, it is fine in steps 3 4 5 to use a
comment just saying "Experimental - leave at default value"
or something to that effect.



ADDING A NEW NON-CONFIGURABLE OPTION

1. Add the option to class CommandLineOnlyOptions.
2. Modify AssemblerOptions::addCommandLineOnlyOptions to reflect the new option,
   making sure to include a default value and at least a minimal help message.
3. Document the option in shasta/docs/CommandLineOptions.html.
4. If the option requires validation add it at the appropriate place in
   shasta/srcMain/main.cpp.

For options not ready for end users, it is fine in steps 3 4 5 to use a
comment just saying "Experimental - leave at default value"
or something to that effect.

*******************************************************************************/

// Boost libraries.
#include <boost/program_options.hpp>

// Standard library.
#include "iostream.hpp"
#include "string.hpp"
#include "vector.hpp"

namespace shasta {
    class AssemblerOptions;
    class AssemblyOptions;
    class CommandLineOnlyOptions;
    class KmersOptions;
    class Mode3AssemblyOptions;
    class ReadsOptions;

}



// Options only allowed on the command line and not in the configuration file.
class shasta::CommandLineOnlyOptions {
public:
    string configName;
    vector <string> inputFileNames;
    string assemblyDirectory;
    string command;
    string memoryMode;
    string memoryBacking;
    uint32_t threadCount;
    string exploreAccess;
    uint16_t port;
};



// Options in the [Reads] section of the configuration file.
// Can also be entered on the command line with option names
// beginning with "Reads.".
class shasta::ReadsOptions {
public:
    int minReadLength;
};



// Options in the [Kmers] section of the configuration file.
// Can also be entered on the command line with option names
// beginning with "Kmers.".
class shasta::KmersOptions {
public:
    int k;
    double probability;
};



// Assembly options that are specific to Mode 3 assembly.
// See source code in the mode3 namespace
// (source files with a mode3- prefix) for more information
class shasta::Mode3AssemblyOptions {
public:

    uint64_t minAnchorCoverage;
    uint64_t maxAnchorCoverage;

    // Options used by class mode3::LocalAssembly
    class LocalAssemblyOptions {
    public:

        // The estimated offset gets extended by this ratio to
        // decide how much to extend reads that only appear in edgeIdA or edgeIdB.
        double estimatedOffsetRatio;

        // Vertex sampling rate, used to set minVertexCoverage.
        // Only used if minVertexCoverage is 0 on input to
        // mode3::LocalAssembly constructor.
        double vertexSamplingRate;

        // Alignment parameters.
        int64_t matchScore;
        int64_t mismatchScore;
        int64_t gapScore;

        // Number of bases (not markers) that can be skipped by an alignment.
        uint64_t maxSkipBases;

        // The maximum tolerated length drift of each read.
        // Used to compute the band for banded alignments.
        double maxDrift;

        // Minimum half band, in markers.
        uint64_t minHalfBand;

        // Minimum ration of scorew to best possible score for
        // an alignment to be used.
        double minScoreRatio;

        // The maximum length of an MSA alignment we are willing to compute.
        uint64_t maxMsaLength;

    };
    LocalAssemblyOptions localAssemblyOptions;
};



// Options in the [Assembly] section of the configuration file.
// Can also be entered on the command line with option names
// beginning with "Assembly.".
class shasta::AssemblyOptions {
public:
    Mode3AssemblyOptions mode3Options;
};



class shasta::AssemblerOptions {
public:

    // Object containing the options.
    CommandLineOnlyOptions commandLineOnlyOptions;
    ReadsOptions readsOptions;
    KmersOptions kmersOptions;
    AssemblyOptions assemblyOptions;

    // Constructor from a command line.
    // If the command line includes a --config option,
    // the specified built-in configuration or configuration file
    // is used to fill the AssemblyOptions,
    // but values specified on the command line take precedence.
    AssemblerOptions(int argumentCount, const char** arguments);

    // Constructor from a configuration file.
    // This only fills in the configurable options specified in
    // the given configuration file. Command line only options
    // are left at their defaults.
    AssemblerOptions(const string& fileName);

    // Add configurable options to the Boost option description object.
    void addCommandLineOnlyOptions();
    void addConfigurableOptions();

    // Boost program_options library objects.
    boost::program_options::options_description commandLineOnlyOptionsDescription;
    boost::program_options::options_description configurableOptionsDescription;
    boost::program_options::options_description allOptionsDescription;

    // This one is the same as allOptionsDescription, with
    // "--invalidOption" added to capture invalid positional options.
    vector<string> invalidPositionalOptions;
    boost::program_options::options_description allOptionsIncludingInvalidDescription;

};

#endif

