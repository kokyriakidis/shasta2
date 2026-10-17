#include "AssemblyGraph.hpp"
#include "findDescendants.hpp"
#include "Options.hpp"
#include "performanceLog.hpp"
#include "ReadFollowing5.hpp"
#include "ReadFollowing6.hpp"
#include "StrandSeparation1.hpp"
#include "SegmentGraph.hpp"
#include "Tangle.hpp"
#include "timestamp.hpp"
using namespace shasta2;
using namespace ReadFollowing5;
using namespace ReadFollowing6;



// This uses ReadFollowing5 code.
bool AssemblyGraph::readFollowingStrandSymmetric5(
    uint64_t tangleId,
    const Tangle& tangle,
    ostream& html,
    const string& debugOutputBaseName)
{
    SHASTA2_ASSERT(not tangle.entrances.empty());
    SHASTA2_ASSERT(not tangle.exits.empty());

    ReadFollowing5::Graph graph(*this, tangleId, tangle, html, debugOutputBaseName);

    if(graph.updateAssemblyGraph1()) {
        return true;
    } else {
        return graph.updateAssemblyGraph2();
    }
}






void AssemblyGraph::separateStrands1(const string& debugOutputBaseName)
{
    // EXPOSE WHEN CODE STABILIZES.
    const uint64_t maxLength = 100000;

    AssemblyGraph& assemblyGraph = *this;

    SHASTA2_ASSERT(countZeroLengthSegments() == 0);

    performanceLog << timestamp << "AssemblyGraph::separateStrands1 begins: " <<
        debugOutputBaseName << endl;

    // Create tangles as connected components of the AssemblyGraph,
    // computed considering only short Segments.
    // The vertices of each tangle are returned sorted by id
    // so we can do binary searches in them.
    vector< vector<vertex_descriptor> > tangles;
    vector<uint64_t> tangleRc;
    createTanglesBySegmentLength(maxLength, tangles, tangleRc);

    // The self-complementary tangles are strand contacts.
    vector< vector<vertex_descriptor> > strandContacts;
    for(uint64_t tangleId=0; tangleId<tangles.size(); tangleId++) {
        if(tangleRc[tangleId] == tangleId) {
            strandContacts.push_back(tangles[tangleId]);
        }
    }

    // Write a csv file that can be imported into Bandage to see
    // the tangles.
    vector<uint64_t> strandContactsRc(strandContacts.size());
    std::iota(strandContactsRc.begin(), strandContactsRc.end(), 0);
    writeTangles(strandContacts, strandContactsRc,
        debugOutputBaseName + "-StrandContacts-Bandage.csv");

    // Process each strand contact separately.
    for(uint64_t strandContactId=0; strandContactId<strandContacts.size(); strandContactId++) {
        const StrandSeparation1::StrandContact strandContact(
            assemblyGraph,
            strandContacts[strandContactId],
            debugOutputBaseName,
            strandContactId);
    }

    removeIsolatedVertices();

    performanceLog << timestamp << "AssemblyGraph::separateStrands1 ends: " <<
        debugOutputBaseName << endl;
}



// This uses ReadFollowing6 code.
bool AssemblyGraph::readFollowingStrandSymmetric6(
    [[maybe_unused]] uint64_t tangleId,
    [[maybe_unused]] const shasta2::Tangle& shasta2Tangle,
    [[maybe_unused]] ostream& html,
    const string& debugOutputBaseName)
{
    SHASTA2_ASSERT(not shasta2Tangle.isSelfComplementary());

    ReadFollowing6::Tangle tangle(
        *this,
        shasta2Tangle.tangleVertices,
        html,
        tangleId,
        debugOutputBaseName);

    return false;
}
