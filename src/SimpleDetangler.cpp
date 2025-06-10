#include "SimpleDetangler.hpp"
#include "Tangle.hpp"
#include "TangleMatrix.hpp"
using namespace shasta;

#include <boost/graph/iteration_macros.hpp>



bool SimpleDetangler::operator()(Tangle& tangle, bool doDetangle)
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

    const uint64_t entranceCount = tangleMatrix.entrances.size();
    const uint64_t exitCount = tangleMatrix.exits.size();

    if(debug) {
        cout << "Working on a tangle with " << entranceCount << " entrances and " << exitCount << " exits." << endl;
        cout << "Entrances:";
        for(const auto& entrance: tangleMatrix.entrances) {
            cout << " " << tangle.assemblyGraph[entrance.e].id;
        }
        cout << endl;
        cout << "Exits:";
        for(const auto& exit: tangleMatrix.exits) {
            cout << " " << tangle.assemblyGraph[exit.e].id;
        }
        cout << endl;
        cout << "Tangle matrix:" << endl;
        for(uint64_t iEntrance=0; iEntrance<entranceCount; iEntrance++) {
            for(uint64_t iExit=0; iExit<exitCount; iExit++) {
                const uint64_t coverage = tangleMatrix.tangleMatrix[iEntrance][iExit].size();
                cout << coverage << " ";
            }
            cout << endl;
        }

        cout << "Internal edges of this tangle:";
        for(const AssemblyGraph::vertex_descriptor v: tangle.tangleVertices) {
            BGL_FORALL_OUTEDGES(v, e, tangle.assemblyGraph, AssemblyGraph) {
                cout << " " << tangle.assemblyGraph[e].id;
            }
        }
        cout << endl;
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
        if(debug) {
            cout << "Not detangling because ambiguous entries are present." << endl;
        }
        return false;
    }

    // If all entrances are connected to all exits, don't do it.
    if(significantEntries.size() == tangleMatrix.entrances.size() * tangleMatrix.exits.size()) {
        if(debug) {
            cout << "Not detangling because all entrances are connected to all exits." << endl;
        }
        return false;
    }

    // If there are no significant entries, don't detangle.
    if(significantEntries.empty()) {
        if(debug) {
            cout << "Not detangling because no significant entries are present." << endl;
        }
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
            if(debug) {
                cout << "Not detangling because not all entrances are connected to at least one exit." << endl;
            }
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
            if(debug) {
                cout << "Not detangling because not all exits are connected to at least one enrance." << endl;
            }
            return false;
        }
    }

    // All good, detangle using the significant entries.
    for(const auto& p: significantEntries) {
        tangle.connect(p.first, p.second);
    }

    if(debug) {
        cout << "Detangling." << endl;
    }

    if(doDetangle) {
        tangle.detangle();
    }

    return true;
}
