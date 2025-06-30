#include "SimpleDetangler.hpp"
#include "Tangle.hpp"
#include "TangleMatrix.hpp"
using namespace shasta;

#include <boost/graph/iteration_macros.hpp>



bool SimpleDetangler::operator()(Tangle& tangle)
{
    const TangleMatrix& tangleMatrix = *(tangle.tangleMatrix);

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

    // If we have ambiguous elements, don't do anything.
    if(not ambiguousEntries.empty()) {
        if(debug) {
            cout << "Not detangling because the tangle matrix has ambiguous elements." << endl;
        }
        return false;
    }

    // If there are no insignificant elements, don't do anything.
    if(insignificantEntries.empty()) {
        if(debug) {
            cout << "Not detangling because the tangle matrix has no insignificant elements." << endl;
        }
        return false;
    }

    // All good, detangle using the significant entries.
    for(const auto& p: significantEntries) {
        tangle.connect(p.first, p.second);
    }

    if(debug) {
        cout << "Detangling." << endl;
    }
    tangle.detangle();

    return true;
}
