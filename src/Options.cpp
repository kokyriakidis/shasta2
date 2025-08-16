// Shasta.
#include "Options.hpp"
using namespace shasta;

// Standard library.
#include <array.hpp>
#include "stdexcept.hpp"
#include <thread>



Options::Options(int argc, char** argv) :
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

    // Adjust the numbers of threads, if necessary.
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }
}



void Options::addOptions()
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

    add_option("--max-marker-error-rate", maxMarkerErrorRate,
        "Maximum marker error rate. Reads with a higher marker error rate are not used."
        )->capture_default_str();

    add_option("--min-anchor-coverage", minAnchorCoverage,
        "Minimum anchor coverage."
        )->capture_default_str();

    add_option("--max-anchor-coverage", maxAnchorCoverage,
        "Maximum anchor coverage."
        )->capture_default_str();

    add_option("--max-anchor-repeat-length", maxAnchorRepeatLength,
        "Maximum number of copies of repeats of period 1, 2, 3,... allowed in an anchor sequence."
        )->capture_default_str();

    add_option("--min-anchor-graph-edge-coverage", minAnchorGraphEdgeCoverage,
        "Minimum anchor graph edge coverage."
        )->capture_default_str();

    add_option("--min-anchor-graph-continue-read-following-count", minAnchorGraphContinueReadFollowingCount,
        "Coverage threshold to continue read following during anchor graph edge creation."
        )->capture_default_str();

    add_option("--transitive-reduction-max-edge-coverage", transitiveReductionMaxEdgeCoverage,
        "Maximum coverage of an AnchorGraph edge subject to removal during transitive reduction."
        )->capture_default_str();

    add_option("--transitive-reduction-max-distance", transitiveReductionMaxDistance,
        "Maximum distance for transitive reduction of the AnchorGraph."
        )->capture_default_str();

    add_option("--a-drift", aDrift,
        "Constant allowed offset for offset consistency."
        )->capture_default_str();

    add_option("--b-drift", bDrift,
        "Constant factor for offset consistency."
        )->capture_default_str();

    add_option("--bubble-cleanup-max-bubble-length",
       bubbleCleanupMaxBubbleLength,
       "Maximum bubble length for bubble cleanup."
       )->capture_default_str();

    add_option("--bubble-cleanup-min-common-count",
       bubbleCleanupMinCommonCount,
       "Minimum number of common oriented reads for bubble cleanup."
       )->capture_default_str();

    add_option("--clustering-min-jaccard",
        clusteringMinJaccard,
       "Minimum Jaccard similarity for oriented read clustering."
       )->capture_default_str();

    add_option("--find-superbubbles-max-distance",
        findSuperbubblesMaxDistance,
       "Maximum distance (number of BFS edges) when finding superbubbles."
       )->capture_default_str();

    add_option("--simplify-superbubble-min-coverage",
        simplifySuperbubbleMinCoverage,
       "Minimum coverage when simplifying superbubbles."
       )->capture_default_str();

    add_option("--simplify-superbubble-max-offset",
        simplifySuperbubbleMaxOffset,
       "Maximum offset when simplifying superbubbles."
       )->capture_default_str();

    add_option("--phasing-distance",
        phasingDistance,
       "Maximum bubble-bubble distance (number of edges) for phasing."
       )->capture_default_str();

    add_option("--phasing-min-degree",
        phasingMinDegree,
       "Minimum degree for phasing graph vertices."
       )->capture_default_str();

    add_option("--phasing-min-coverage",
        phasingMinCoverage,
       "Minimum coverage generated by phasing."
       )->capture_default_str();

    add_option("--detangle-min-common-coverage",
        detangleMinCommonCoverage,
        "Minimum common coverage for detangling."
        )->capture_default_str();

    add_option("--detangle-low-coverage-threshold",
        detangleLowCoverageThreshold,
        "Low coverage threshold for detangling."
        )->capture_default_str();

    add_option("--detangle-high-coverage-threshold",
        detangleHighCoverageThreshold,
        "High coverage threshold for detangling."
        )->capture_default_str();

    add_option("--detangle-max-trim",
        detangleMaxTrim,
        "Maximum trim for detangling (number of assembly steps)."
        )->capture_default_str();


    add_option("--detangle-epsilon",
        detangleEpsilon,
        "Epsilon value for chi-square detangling."
        )->capture_default_str();

    add_option("--detangle-maxLogP",
        detangleMaxLogP,
        "MaxLogP value for chi-square detangling."
        )->capture_default_str();

    add_option("--detangle-minLogPDelta",
        detangleMinLogPDelta,
        "MinLogPDelta value for chi-square detangling."
        )->capture_default_str();

    add_option("--prune-length",
        pruneLength,
        "Maximum length of a hanging segments to be pruned."
        )->capture_default_str();

    add_option("--prune-iteration-count",
        pruneIterationCount,
        "Maximum number of pruning iterations."
        )->capture_default_str();

    add_option("--max-abpoa-length", maxAbpoaLength,
        "Maximum MSA length for abpoa (switch to poasta above that)."
        )->capture_default_str();
}



// Constructor from a configuration file.
Options::Options(const string& fileName)
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

    // Adjust the numbers of threads, if necessary.
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }
}



void Options::write(ostream& s) const
{
    s << config_to_str(true,true);
}
