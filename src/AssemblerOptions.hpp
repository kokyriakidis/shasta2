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

    // Function to convert a bool to True or False for better
    // compatibility with Python scripts.
    string convertBoolToPythonString(bool);
}



// Options only allowed on the command line and not in the configuration file.
class shasta::CommandLineOnlyOptions {
public:
    string configName;
    vector <string> inputFileNames;
    vector <string> anchorFileNames;
    string assemblyDirectory;
    string command;
    string memoryMode;
    string memoryBacking;
    uint32_t threadCount;
    bool suppressStdoutLog;
    string exploreAccess;
    uint16_t port;
    string alignmentsPafFile;
    bool saveBinaryData;
};



// Options in the [Reads] section of the configuration file.
// Can also be entered on the command line with option names
// beginning with "Reads.".
class shasta::ReadsOptions {
public:
    uint64_t representation;    // 0 = Raw, 1=RLE
    int minReadLength;
    bool noCache;

    // String to control handling of duplicate reads.
    // Can be one of:
    // useAllCopies
    // useOneCopy
    // useNone
    // forbid
    // See ReadFlags.hpp for the meaning of each option.
    string handleDuplicates;

    void write(ostream&) const;
};



// Options in the [Kmers] section of the configuration file.
// Can also be entered on the command line with option names
// beginning with "Kmers.".
class shasta::KmersOptions {
public:
    int generationMethod;
    int k;
    double probability;
    double enrichmentThreshold;
    uint64_t distanceThreshold;
    string file;
    string globalFrequencyOverrideDirectory;
    void write(ostream&) const;
};



// Assembly options that are specific to Mode 3 assembly.
// See source code in the mode3 namespace
// (source files with a mode3- prefix) for more information
class shasta::Mode3AssemblyOptions {
public:

    string anchorCreationMethod;

    uint64_t minAnchorCoverage;
    uint64_t maxAnchorCoverage;
    double minAnchorCoverageMultiplier;
    double maxAnchorCoverageMultiplier;

    // Options used to clean up the PrimaryGraph.
    class PrimaryGraphOptions {
    public:

        // Parameter to control removal of weak edges.
        double maxLoss;

        // Parameters to control removal of cross edges.
        uint64_t crossEdgesLowCoverageThreshold;
        uint64_t crossEdgesHighCoverageThreshold;

        void write(ostream&) const;
    };
    PrimaryGraphOptions primaryGraphOptions;



    class AssemblyGraphOptions {
    public:

        // Detangle tolerances.
        uint64_t detangleToleranceLow;
        uint64_t detangleToleranceHigh;

        // Bayesian model.
        double epsilon;
        double minLogP;

        // Other thresholds used by the mode3::AssemblyGraph
        uint64_t longBubbleThreshold;
        double phaseErrorThreshold;
        double bubbleErrorThreshold;
        uint64_t bubbleCleanupMaxOffset;
        uint64_t chainTerminalCommonThreshold;
        uint64_t superbubbleLengthThreshold1;
        uint64_t superbubbleLengthThreshold2;
        uint64_t superbubbleLengthThreshold3;
        uint64_t superbubbleLengthThreshold4;
        uint64_t pruneLength;

        void write(ostream&) const;
    };
    AssemblyGraphOptions assemblyGraphOptions;



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

        void write(ostream&) const;
    };
    LocalAssemblyOptions localAssemblyOptions;

    void write(ostream&) const;
};



// Options in the [Assembly] section of the configuration file.
// Can also be entered on the command line with option names
// beginning with "Assembly.".
class shasta::AssemblyOptions {
public:
    uint64_t mode;
    int crossEdgeCoverageThreshold;
    int markerGraphEdgeLengthThresholdForConsensus;
    string consensusCallerString;
    string consensusCaller;
    bool storeCoverageData;
    int storeCoverageDataCsvLengthThreshold;
    bool writeReadsByAssembledSegment;
    uint64_t pruneLength;

    // Options that control detangling.
    int detangleMethod;
    uint64_t detangleDiagonalReadCountMin;
    uint64_t detangleOffDiagonalReadCountMax;
    double detangleOffDiagonalRatio;

    // Options that control iterative assembly.
    bool iterative;
    uint64_t iterativeIterationCount;
    int64_t iterativePseudoPathAlignMatchScore;
    int64_t iterativePseudoPathAlignMismatchScore;
    int64_t iterativePseudoPathAlignGapScore;
    double iterativeMismatchSquareFactor;
    double iterativeMinScore;
    uint64_t iterativeMaxAlignmentCount;
    uint64_t iterativeBridgeRemovalIterationCount;
    uint64_t iterativeBridgeRemovalMaxDistance;

    // Mode 3 assembly options.
    Mode3AssemblyOptions mode3Options;

    void write(ostream&) const;

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

    // Write the options as a config file.
    void write(ostream&) const;

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

