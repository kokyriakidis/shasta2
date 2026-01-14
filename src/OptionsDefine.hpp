// General options.
SHASTA2_OPTION_DEFINE(
    string, command, "--command", "assemble",
    "Shasta2 command.")

SHASTA2_OPTION_DEFINE(
    string, assemblyDirectory, "--output", "ShastaRun",
    "Assembly directory for output.")

SHASTA2_OPTION_DEFINE(
    string, memoryMode, "--memory-mode", "anonymous",
    "Memory mode. Specify anonymous or filesystem.")

SHASTA2_OPTION_DEFINE(
    string, memoryBacking, "--memory-backing", "4K",
    "Memory backing. Specify disk, 4K, or 2M.")

SHASTA2_OPTION_DEFINE(
    uint64_t, threadCount, "--threads", 0,
    "Number of threads (0 to use hardware_concurrency).")

SHASTA2_OPTION_DEFINE(
    string, externalAnchorsName, "--external-anchors-name", "",
    "External anchors name. Must be an absolute path without the .toc/.data extensions.")



// Http server options.
SHASTA2_OPTION_DEFINE(
    string, exploreAccess, "--explore-access", "user",
        "Specifies accessibility of Shasta2 http server.\n"
        "Specify user, local, or unrestricted.")

SHASTA2_OPTION_DEFINE(
    uint16_t, port, "--port", 17100,
    "Port number for http server.")



// Reads.
SHASTA2_OPTION_DEFINE(
    uint64_t, minReadLength, "--min-read-length", 0,
    "Read length cutoff.")



// Markers.
SHASTA2_OPTION_DEFINE(
    uint64_t, k, "--k", 60,
    "Marker length.")

SHASTA2_OPTION_DEFINE(
    double, markerDensity, "--marker-density", 0.05,
    "Marker density.")

SHASTA2_OPTION_DEFINE(
    double, maxMarkerErrorRate, "--max-marker-error-rate", 0.5,
    "Maximum marker error rate. Reads with a higher marker error rate are not used.")



// Anchors.
SHASTA2_OPTION_DEFINE(
    uint64_t, minAnchorCoverage, "--min-anchor-coverage", 10,
    "Minimum anchor coverage.")

SHASTA2_OPTION_DEFINE(
    uint64_t, maxAnchorCoverage, "--max-anchor-coverage", 60,
    "Maximum anchor coverage.")



// Options that control detangling.
SHASTA2_OPTION_DEFINE(
    double, detangleEpsilon, "--detangle-epsilon", 0.05,
    "Epsilon value for likelihood ratio detangling.")

SHASTA2_OPTION_DEFINE(
    double, detangleMaxLogP, "--detangle-maxLogP", 30.,
    "MaxLogP value for likelihood ratio detangling.")

SHASTA2_OPTION_DEFINE(
    double, detangleMinLogPDelta, "--detangle-minLogPDelta", 10.,
    "MinLogPDelta value for likelihood ratio detangling.")

SHASTA2_OPTION_DEFINE(
    uint64_t ,representativeRegionStepCount, "--representative-region-step-count", 10,
    "Number of steps for the representative region of a segment "
    "used in detangling and read following.")



// Options that control read following.
SHASTA2_OPTION_DEFINE(
	uint64_t, readFollowingMinCommonCount, "--read-following-min-common-count", 6,
	"Minimum number of common oriented reads for read following.")

SHASTA2_OPTION_DEFINE(
	double, readFollowingMinCorrectedJaccard, "--read-following-min-corrected-jaccard", 0.7,
	"Minimum corrected Jaccard for read following.")

SHASTA2_OPTION_DEFINE(
	uint32_t, readFollowingPruneLength, "--read-following-prune-length", 100000,
	"Prune length for read following.")

SHASTA2_OPTION_DEFINE(
	uint64_t, readFollowingSegmentLengthThreshold, "--read-following-segment-length-threshold", 500000,
	"Segment length threshold for read following.")

