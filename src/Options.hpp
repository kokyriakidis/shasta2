#pragma once

// Shasta.
#include <cstdint.hpp>
#include <iosfwd.hpp>
#include <string.hpp>
#include <vector.hpp>

// CLI11. Requires package libcli11-dev.
#include "CLI/CLI.hpp"

namespace shasta {
    class Options;
}


class shasta::Options : public CLI::App {
public:

    // Constructor from command line options.
    Options(int argumentCount, char** arguments);

    // Constructor from a configuration file.
    Options(const string& fileName);

    // Write a configuration file.
    void write(ostream&) const;

    bool isHelp = false;

    string configName;
    vector <string> inputFileNames;
    string assemblyDirectory = "ShastaRun";
    string command = "assemble";
    string memoryMode = "anonymous";
    string memoryBacking= "4K";
    uint64_t threadCount = 0;
    string exploreAccess = "user";
    uint16_t port = 17100;

    uint64_t minReadLength = 0;

    uint64_t k = 60;
    double markerDensity = 0.05;

    double maxMarkerErrorRate = 0.5;

    string externalAnchorsName;
    uint64_t minAnchorCoverage = 10;
    uint64_t maxAnchorCoverage = 60;

    // An anchor is not generated if its sequence contains an exact repeat
    // consisting of n copies of a unit of length (period) p, if
    // n > maxAnchorRepeatLength[p-1].
    // So for example:
    // - maxAnchorRepeatLength[0] is the maximum allowed length of
    //   a homopolymer run.
    // - maxAnchorRepeatLength[1] is the maximum allowed number of
    //   copies of a repeat with period 2 (e. g. ATATAT).
    //   Note this is the number of copies, not number of bases.
    //   So if maxAnchorRepeatLength[1] is 3, the anchor is not
    //   generated if it contains a 2-base repeat with more than 3 copies
    //   (a total 6 bases).
    vector<uint64_t> maxAnchorRepeatLength = {6, 4, 4, 4, 4};

    // Options controlling creation of the AnchorGraph.
    uint64_t minAnchorGraphEdgeCoverage = 6;
    uint64_t minAnchorGraphContinueReadFollowingCount = 10;
    uint64_t transitiveReductionMaxEdgeCoverage = 10;
    uint64_t transitiveReductionMaxDistance = 10;

    double aDrift = 300.;
    double bDrift = 0.01;

    uint64_t simplifyMaxIterationCount = 3;

    uint64_t bubbleCleanupMaxBubbleLength = 10000;
    uint64_t bubbleCleanupMinCommonCount = 6;
    double clusteringMinJaccard = 0.7;

    uint64_t findSuperbubblesMaxDistance = 10;
    uint64_t simplifySuperbubbleMinCoverage = 4;
    uint64_t simplifySuperbubbleMaxOffset = 30000;

    // Options that control phasing.
    uint64_t phasingDistance = 12;
    uint64_t phasingMinDegree = 2;
    uint64_t phasingMinCoverage = 4;

    // Options that control detangling.
    double detangleEpsilon = 0.05;
    double detangleMaxLogP = 30.;
    double detangleMinLogPDelta = 10.;
    uint64_t detangleMaxIterationCount = 10;
    uint64_t detangleMaxCrossEdgeLength = 10000;
    uint64_t detangleMinCoverage = 3;
    uint64_t representativeRegionStepCount = 10;

    // Options that control pruning.
    uint64_t pruneLength = 50000;
    uint64_t pruneIterationCount = 3;

	// Options defined in OptionsDefine.hpp
	#define SHASTA2_OPTION_DEFINE(type, name, optionName, defaultValue, description) \
		type name = defaultValue;
	#include "OptionsDefine.hpp"
	#undef SHASTA2_OPTION_DEFINE

private:
    void addOptions();
};


