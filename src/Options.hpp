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

    uint64_t minAnchorCoverage = 10;
    uint64_t maxAnchorCoverage = 60;
    uint64_t maxAnchorHomopolymerLength = 10;
    uint64_t minAnchorGraphEdgeCoverage = 6;
    uint64_t minAnchorGraphContinueReadFollowingCount = 10;

    // Transitive reduction of the AnchorGraph.
    uint64_t transitiveReductionMaxEdgeCoverage = 10;
    uint64_t transitiveReductionMaxDistance = 10;

    double aDrift = 300.;
    double bDrift = 0.01;

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
    uint64_t detangleMinCommonCoverage = 3;
    uint64_t detangleLowCoverageThreshold = 1;
    uint64_t detangleHighCoverageThreshold = 4;
    uint64_t detangleMaxTrim = 10;
    double detangleEpsilon = 0.05;
    double detangleMaxLogP = 30.;
    double detangleMinLogPDelta = 10.;

    // Options that control pruning.
    uint64_t pruneLength = 50000;
    uint64_t pruneIterationCount = 3;

    // Maximum MSA length for abpoa (switch to poasta above that).
    uint64_t maxAbpoaLength = 5000;

private:
    void addOptions();
};


