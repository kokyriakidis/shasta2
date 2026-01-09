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
    uint64_t, detangleMaxIterationCount, "--detangle-max-iteration-count", 10,
    "Maximum number of detangling iterations at each iteration of assembly graph simplify.")

SHASTA2_OPTION_DEFINE(
    uint64_t, detangleMaxCrossEdgeLength, "--detangle-max-cross-edge-length", 10000,
    "Maximum cross-edge length for detangling.")

SHASTA2_OPTION_DEFINE(
    uint64_t, detangleMinCoverage, "--detangle-min-coverage", 3,
    "Minimum coverage that phasing and detangling are allowed to generate.")

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



SHASTA2_OPTION_DEFINE(
	uint64_t, maxAbpoaLength, "--max-abpoa-length", 5000,
	"Maximum MSA length for abpoa (switch to poasta above that).")
