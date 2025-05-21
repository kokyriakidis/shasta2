#include "SimpleDetangler.hpp"
#include "Tangle.hpp"
#include "Tangle2.hpp"
#include "Tangle3.hpp"
#include "TangleMatrix2.hpp"
#include "TangleMatrix3.hpp"
using namespace shasta;



bool SimpleDetangler::operator()(Tangle&)
{
    return false;
}



bool SimpleDetangler::operator()(Tangle2&)
{
    return false;

#if 0
    const TangleMatrix2& tangleMatrix = *(tangle.tangleMatrix);

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


    if(tangleMatrix.entrances.size() < 2) {
        return false;
    }
    if(tangleMatrix.exits.size() < 2) {
        return false;
    }

    // Check that each entrance will get a connection with at least one exit.
    for(uint64_t iEntrance=0; iEntrance<tangleMatrix.entrances.size(); iEntrance++) {
        bool isGood = false;
        for(uint64_t iExit=0; iExit<tangleMatrix.exits.size(); iExit++) {
            if(tangleMatrix.tangleMatrix[iEntrance][iExit].size() >= minDetangleCoverage) {
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
        for(uint64_t iEntrance=0; iEntrance<tangleMatrix.entrances.size(); iEntrance++) {
            if(tangleMatrix.tangleMatrix[iEntrance][iExit].size() >= minDetangleCoverage) {
                isGood = true;
                break;
            }
        }
        if(not isGood) {
            return false;
        }
    }


    // Do the detangling.
    // Connect the pairs with minDetangleCoverage.
    for(uint64_t iEntrance=0; iEntrance<tangleMatrix.entrances.size(); iEntrance++) {
        for(uint64_t iExit=0; iExit<tangleMatrix.exits.size(); iExit++) {
            if(tangleMatrix.tangleMatrix[iEntrance][iExit].size() >= minDetangleCoverage) {
                tangle.connect(iEntrance, iExit);
            }
        }

    }
    tangle.detangle();

    return true;
#endif
}



bool SimpleDetangler::operator()(Tangle3& tangle)
{
    const TangleMatrix3& tangleMatrix = *(tangle.tangleMatrix);

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

    // For now only do the 2x2 case.
    if(tangleMatrix.entrances.size() < 2) {
        return false;
    }
    if(tangleMatrix.exits.size() < 2) {
        return false;
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


    // We can generate a connection for each significant entry.

    // If all entrances are connected to all exits, don't do it.
    if(significantEntries.size() == tangleMatrix.entrances.size() * tangleMatrix.exits.size()) {
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

    // All good, detangle using the significant entries.
    for(const auto& p: significantEntries) {
        tangle.connect(p.first, p.second);
    }
    tangle.detangle();

    return true;
}
