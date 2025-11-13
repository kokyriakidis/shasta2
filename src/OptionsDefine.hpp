
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
