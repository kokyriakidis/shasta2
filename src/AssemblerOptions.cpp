// Shasta.
#include "AssemblerOptions.hpp"
using namespace shasta;

// Standard library.
#include <array.hpp>
#include "stdexcept.hpp"



shasta::AssemblerOptions::AssemblerOptions(int argc, char** argv) :
    CLI::App("Shasta2. Under development.")
{
    allow_config_extras(false);
    set_config("--config", "", "Configuration file.");

    get_formatter()->column_width(20);

    addOptions();

    try {
        parse(argc, argv);
    } catch(const CLI::ParseError& e) {
         exit(e);
         if(e.get_name() == "CallForHelp") {
             isHelp = true;
         } else {
             throw runtime_error("Error parsing options.");
         }
    }
}



void AssemblerOptions::addOptions()
{
    add_option("--input", inputFileNames,
        "Input fasta or fastq files (uncompressed)."
        )->capture_default_str();

    add_option("--output", assemblyDirectory,
        "Assembly directory for output. Must not exist."
        )->capture_default_str();

    add_option("--command", command,
        "Shasta2 command.")
    ->capture_default_str();

    add_option("--memory-mode", memoryMode,
        "Memory mode. Specify disk, 4K, or 2M.")
    ->capture_default_str();

    add_option("--memory-backing", memoryBacking,
        "Memory backing.Specify anonymous or filesystem."
        )->capture_default_str();

    add_option("--threads", threadCount,
        "Number of threads (0 to use hardware_concurrency)."
        )->capture_default_str();

    add_option("--explore-access", exploreAccess,
        "Specifies accessibility of Shasta2 http server.\n"
        "Specify user, local, or unrestricted."
        )->capture_default_str();

    add_option("--port", port,
        "Port number for http server."
        )->capture_default_str();

    add_option("--min-read-length", minReadLength,
        "Read length cutoff."
        )->capture_default_str();

    add_option("--k", k,
        "Marker length"
        )->capture_default_str();

    add_option("--marker-density", markerDensity,
        "Marker density."
        )->capture_default_str();

    add_option("--min-anchor-coverage", minAnchorCoverage,
        "Minimum anchor coverage."
        )->capture_default_str();

    add_option("--max-anchor-coverage", maxAnchorCoverage,
        "Maximum anchor coverage."
        )->capture_default_str();

    add_option("--max-anchor-homopolymer-length", maxAnchorHomopolymerLength,
        "Maximum homopolymer length allowed on an anchor."
        )->capture_default_str();

    add_option("--min-anchor-graph-edge-coverage-near", minAnchorGraphEdgeCoverageNear,
        "Minimum anchor graph edge coverage for small offset."
        )->capture_default_str();

    add_option("--min-anchor-graph-edge-coverage-far", minAnchorGraphEdgeCoverageFar,
        "Minimum anchor graph edge coverage for large offset."
        )->capture_default_str();

    add_option("--anchor-graph-edge-coverage-fraction-threshold", anchorGraphEdgeCoverageFractionThreshold,
        "Anchor graph edge coverage fraction threshold."
        )->capture_default_str();

    add_option("--a-drift", aDrift,
        "Constant allowed offset for offset consistency."
        )->capture_default_str();

    add_option("--b-drift", bDrift,
        "Constant factor for offset consistency."
        )->capture_default_str();

     add_option("--bubble-cleanup-min-common-count",
        bubbleCleanupMinCommonCount,
        "Minimum number of common oriented reads for bubble cleanup."
        )->capture_default_str();

    add_option("--detangle-min-common-coverage",
        assemblyGraphOptions.detangleMinCommonCoverage,
        "Minimum common coverage for detangling."
        )->capture_default_str();

    add_option("--detangle-low-coverage-threshold",
        assemblyGraphOptions.detangleLowCoverageThreshold,
        "Low coverage threshold for detangling."
        )->capture_default_str();

    add_option("--detangle-high-coverage-threshold",
        assemblyGraphOptions.detangleHighCoverageThreshold,
        "High coverage threshold for detangling."
        )->capture_default_str();

    add_option("--detangle-epsilon",
        assemblyGraphOptions.detangleEpsilon,
        "Epsilon value for chi-square detangling."
        )->capture_default_str();

    add_option("--detangle-maxLogP",
        assemblyGraphOptions.detangleMaxLogP,
        "MaxLogP value for chi-square detangling."
        )->capture_default_str();

    add_option("--detangle-minLogPDelta",
        assemblyGraphOptions.detangleMinLogPDelta,
        "MinLogPDelta value for chi-square detangling."
        )->capture_default_str();

    add_option("--prune-length",
        pruneLength,
        "Maximum length of a hanging segments to be pruned."
        )->capture_default_str();

    add_option("--local-assembly-estimated-offset-ratio",
        localAssemblyOptions.estimatedOffsetRatio,
        "Estimated offset ratio for local assembly."
        )->capture_default_str();

    add_option("--local-assembly-vertex-sampling-rate", localAssemblyOptions.vertexSamplingRate,
        "Vertex sampling rate for local assembly."
        )->capture_default_str();

    add_option("--local-assembly-match-score", localAssemblyOptions.matchScore,
        "Match score for local assembly."
        )->capture_default_str();

    add_option("--local-assembly-mismatch-score", localAssemblyOptions.mismatchScore,
        "Mismatch score for local assembly."
        )->capture_default_str();

    add_option("--local-assembly-gap-score", localAssemblyOptions.gapScore,
        "Gap score for local assembly."
        )->capture_default_str();

    add_option("--local-assembly-max-skip-bases", localAssemblyOptions.maxSkipBases,
        "Maximum skip in bases for local assembly."
        )->capture_default_str();

    add_option("--local-assembly-max-drift-ratio", localAssemblyOptions.maxDrift,
        "Maximum drift ratio for local assembly."
        )->capture_default_str();

    add_option("--local-assembly-min-half-band", localAssemblyOptions.minHalfBand,
        "Minimum half band, in markers, for local assembly."
        )->capture_default_str();

    add_option("--local-assembly-min-score-ratio", localAssemblyOptions.minScoreRatio,
        "Score ratio threshold for discarding alignments\n"
        "for local assembly."
        )->capture_default_str();

    add_option("--local-assembly-max-msa-length", localAssemblyOptions.maxMsaLength,
        "Maximum length of a multiple sequence alignment\n"
        "for local assembly."
        )->capture_default_str();

    add_option("--local-assembly-max-abpoa-length", localAssemblyOptions.maxAbpoaLength,
        "Maximum MSA length for abpoa (switch to poasta above that)."
        )->capture_default_str();
}



// Constructor from a configuration file.
AssemblerOptions::AssemblerOptions(const string& fileName)
{

    allow_config_extras(false);
    set_config("--config", "", "Configuration file.");
    get_formatter()->column_width(20);

    addOptions();

    // Construct arguments "--config fileName".
    int argc = 3;
    const string name = "shasta2";
    const string keyword = "--config";
    char* arg0 = const_cast<char*>(name.c_str());
    char* arg1 = const_cast<char*>(keyword.c_str());
    char* arg2 = const_cast<char*>(fileName.c_str());
    const array<char*, 3> argv = {arg0, arg1, arg2};

    try {
        parse(argc, &argv.front());
    } catch(const CLI::ParseError& e) {
         exit(e);
         throw runtime_error("Error parsing options.");
    }
}



void AssemblerOptions::write(ostream& s) const
{
    s << config_to_str(true,true);
}
