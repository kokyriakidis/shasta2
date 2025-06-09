#include "SimpleDetangler.hpp"
#include "Tangle.hpp"
#include "TangleMatrix.hpp"
using namespace shasta;



bool SimpleDetangler::operator()(Tangle& tangle)
{
    const TangleMatrix& tangleMatrix = *(tangle.tangleMatrix);

    // Check common coverage on all entrances and exits.
    for(const auto& entrance: tangleMatrix.entrances) {
        if(entrance.commonCoverage < minCommonCoverage) {
            return false;
        }
    }
    for(const auto& exit: tangleMatrix.exits) {
        if(exit.commonCoverage < minCommonCoverage) {
            return false;
        }
    }

    // Gather entries by type.
    vector< pair<uint64_t, uint64_t> > insignificantEntries;
    vector< pair<uint64_t, uint64_t> > ambiguousEntries;
    vector< pair<uint64_t, uint64_t> > significantEntries;
    for(uint64_t iEntrance=0; iEntrance<tangleMatrix.entrances.size(); iEntrance++) {
        for(uint64_t iExit=0; iExit<tangleMatrix.exits.size(); iExit++) {
            const uint64_t coverage = tangleMatrix.tangleMatrix[iEntrance][iExit].size();

            if(coverage <= detangleLowThreshold) {
                insignificantEntries.emplace_back(iEntrance, iExit);
            } else if(coverage >= detangleHighThreshold) {
                significantEntries.emplace_back(iEntrance, iExit);
            } else {
                ambiguousEntries.emplace_back(iEntrance, iExit);
            }
        }
    }

    // If we have ambiguous entries, give up.
    if(not ambiguousEntries.empty()) {
        return false;
    }

    // If all entrances are connected to all exits, don't do it.
    if(significantEntries.size() == tangleMatrix.entrances.size() * tangleMatrix.exits.size()) {
        return false;
    }

    // If there are no significant entries, don't detangle.
    if(significantEntries.empty()) {
        return false;
    }

    // Check that each entrance will get a connection with at least one exit.
    for(uint64_t iEntrance=0; iEntrance<tangleMatrix.entrances.size(); iEntrance++) {
        bool isGood = false;
        for(const auto& p: significantEntries) {
            if(p.first == iEntrance) {
                isGood = true;
                break;
            }
        }
        if(not isGood) {
            return false;
        }
    }

    // Check that each exit will get a connection with at least one entrance.
    for(uint64_t iExit=0; iExit<tangleMatrix.exits.size(); iExit++) {
        bool isGood = false;
        for(const auto& p: significantEntries) {
            if(p.second == iExit) {
                isGood = true;
                break;
            }
        }
        if(not isGood) {
            return false;
        }
    }

#if 0
    // Check the base offsets of the significant entries. If any of them are too long,
    // don't detangle.
    for(const auto& p: significantEntries) {
        const AnchorPair& anchorPair = tangleMatrix.tangleMatrix[p.first][p.second];
        const uint64_t offset = anchorPair.getAverageOffset(tangle.assemblyGraph.anchors);
        if(debug) {
            cout << "Base offset for " << p.first << " " << p.second << " is " << offset << endl;
        }
        if(offset > maxBaseOffset) {
            return false;
        }
    }
#endif

    // All good, detangle using the significant entries.
    for(const auto& p: significantEntries) {
        tangle.connect(p.first, p.second);
    }
    tangle.detangle();

    return true;
}
