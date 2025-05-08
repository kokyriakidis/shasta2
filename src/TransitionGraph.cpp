// Shasta.
#include "TransitionGraph.hpp"
#include "Anchor.hpp"
#include "deduplicate.hpp"
#include "Markers.hpp"
#include "orderPairs.hpp"
#include "performanceLog.hpp"
#include "timestamp.hpp"
using namespace shasta;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>

// Standard library.
#include <fstream.hpp>
#include <tuple.hpp>




TransitionGraph::TransitionGraph(
    const Anchors& anchors,
    const AnchorGraph& anchorGraph,
    uint64_t minTransitionGraphEdgeCoverage)
{
    TransitionGraph& transitionGraph = *this;

    // Create the TransitionGraph vertices, one for each AnchorGraph edge.
    uint64_t nextVertexId = 0;;
    BGL_FORALL_EDGES(eA, anchorGraph, AnchorGraph) {
        add_vertex(TransitionGraphVertex(anchorGraph[eA], nextVertexId++), transitionGraph);
    }

    // Compute journeys.
    computeJourneys(anchors);

    // Gather pairs of consecutive vertices in the journeys.
    performanceLog << timestamp << "Finding journey pairs begins." << endl;
    class Info {
    public:
        vertex_descriptor v0;
        vertex_descriptor v1;
        OrientedReadId orientedReadId;

        Info(
            vertex_descriptor v0,
            vertex_descriptor v1,
            OrientedReadId orientedReadId) :
            v0(v0), v1(v1), orientedReadId(orientedReadId) {}

        bool operator<(const Info& that) const
        {
            return tie(v0, v1, orientedReadId) < tie(that.v0, that.v1, that.orientedReadId);
        }
    };
    vector<Info> infos;
    for(uint64_t j=0; j<journeys.size(); j++) {
        const OrientedReadId orientedReadId = OrientedReadId::fromValue(ReadId(j));
        const vector<vertex_descriptor>& journey = journeys[j];
        for(uint64_t i1=1; i1<journey.size(); i1++) {
            const uint64_t i0 = i1 - 1;
            const vertex_descriptor v0 = journey[i0];
            const vertex_descriptor v1 = journey[i1];
            infos.emplace_back(v0, v1, orientedReadId);
        }
    }
    performanceLog << timestamp << "Finding journey pairs ends." << endl;

    performanceLog << timestamp << "Sort journey pairs begins." << endl;
    sort(infos.begin(), infos.end());
    performanceLog << timestamp << "Sort journey pairs ends." << endl;



    // Loop over streaks of Infos with the same eA0 and eA1.
    // Streaks of length at least minTransitionGraphEdgeCoverage generate
    // a TransitionGraph edge.
    for(uint64_t streakBegin=0; streakBegin!=infos.size(); /* Increment later */) {
        const vertex_descriptor v0 = infos[streakBegin].v0;
        const vertex_descriptor v1 = infos[streakBegin].v1;

        uint64_t streakEnd = streakBegin + 1;
        for(; streakEnd!=infos.size(); streakEnd++) {
            if(infos[streakEnd].v0 != v0) {
                break;
            }
            if(infos[streakEnd].v1 != v1) {
                break;
            }
        }

        // If the streak is sufficiently long, one of two things happen:
        // - If the source and target AnchorPairs of these two
        //   TransitionGraph vertices are adjacent, just generate an edge.
        // - Otherwise, generate a new TransitionGraph vertex to "bridge"
        //   the AnchorPairs, plus two edges.
        if(streakEnd - streakBegin >= minTransitionGraphEdgeCoverage) {

            const TransitionGraphVertex& vertex0 = transitionGraph[v0];
            const TransitionGraphVertex& vertex1 = transitionGraph[v1];

            const AnchorPair& anchorPair0 = vertex0.anchorPair;
            const AnchorPair& anchorPair1 = vertex1.anchorPair;

            if(anchorPair0.anchorIdB == anchorPair1.anchorIdA) {
                add_edge(v0, v1, transitionGraph);
            } else {

                // Create a new TransitionGraph vertex.
                vertex_descriptor v2 = add_vertex(transitionGraph);
                TransitionGraphVertex& vertex2 = transitionGraph[v2];
                vertex2.id = nextVertexId++;
                vertex2.anchorPair.anchorIdA = vertex0.anchorPair.anchorIdB;
                vertex2.anchorPair.anchorIdB = vertex1.anchorPair.anchorIdA;
                for(uint64_t i=streakBegin; i!=streakEnd; i++) {
                    vertex2.anchorPair.orientedReadIds.push_back(infos[i].orientedReadId);
                }

                /*
                cout << "Added bridge vertex to the TransitionGraph between (" <<
                    anchorIdToString(vertex0.anchorPair.anchorIdA) << "," <<
                    anchorIdToString(vertex0.anchorPair.anchorIdB) << ") and (" <<
                    anchorIdToString(vertex1.anchorPair.anchorIdA) << "," <<
                    anchorIdToString(vertex1.anchorPair.anchorIdB) << "), coverage " <<
                    vertex2.anchorPair.orientedReadIds.size() << ", offset " <<
                    vertex2.anchorPair.getAverageOffset(anchors) << endl;
                */


                add_edge(v0, v2, transitionGraph);
                add_edge(v2, v1, transitionGraph);
            }
        }

        // Prepare to process the next streak.
        streakBegin = streakEnd;
    }

    cout << "The TransitionGraph has " << num_vertices(transitionGraph) <<
        " vertices and " << num_edges(transitionGraph) << " edges." << endl;

    // Check that all edges are between adjacent AnchorPairs.
    check();
}



void TransitionGraph::writeGraphviz(const string& fileName) const
{
    ofstream file(fileName);
    writeGraphviz(file);
}



void TransitionGraph::writeGraphviz(ostream& dot) const
{
    const TransitionGraph& transitionGraph = *this;

    dot << "digraph TransitionGraph {\n";

    BGL_FORALL_VERTICES(v, transitionGraph, TransitionGraph) {
        dot << transitionGraph[v].id << ";\n";
    }

    BGL_FORALL_EDGES(e, transitionGraph, TransitionGraph) {
        const vertex_descriptor v0 = source(e, transitionGraph);
        const vertex_descriptor v1 = target(e, transitionGraph);

        dot << transitionGraph[v0].id << "->";
        dot << transitionGraph[v1].id << ";\n";
    }

    dot << "}\n";
}



void TransitionGraph::computeJourneys(const Anchors& anchors)
{
    const TransitionGraph& transitionGraph = *this;
    const uint64_t orientedReadCount = anchors.markers.size();

    // Construct a vector of pairs (ordinal, vertex_descriptor) for each OrientedReadId;
    vector< vector< pair<uint32_t, vertex_descriptor> > > ordinalTable(orientedReadCount);
    vector< pair<uint32_t, uint32_t> > ordinals;
    BGL_FORALL_VERTICES(v, transitionGraph, TransitionGraph) {
        const AnchorPair& anchorPair = transitionGraph[v].anchorPair;
        anchorPair.getOrdinals(anchors, ordinals);
        SHASTA_ASSERT(ordinals.size() == anchorPair.orientedReadIds.size());
        for(uint64_t i=0; i<anchorPair.orientedReadIds.size(); i++) {
            const OrientedReadId orientedReadId = anchorPair.orientedReadIds[i];
            const uint32_t ordinal0 = ordinals[i].first;
            ordinalTable[orientedReadId.getValue()].push_back(make_pair(ordinal0, v));
        }
    }

    // Sort them by ordinals.
    for(vector< pair<uint32_t, vertex_descriptor> >& u: ordinalTable) {
        sort(u.begin(), u.end(), OrderPairsByFirstOnly<uint32_t, vertex_descriptor>());
    }

    // Now we can create the  journeys.
    journeys.clear();
    journeys.resize(orientedReadCount);
    for(uint64_t i=0; i<orientedReadCount; i++) {
        vector<vertex_descriptor>& journey = journeys[i];
        const vector< pair<uint32_t, vertex_descriptor> >& u = ordinalTable[i];
        for(const auto& p: u) {
            journey.push_back(p.second);
        }
    }

}


// This check that the AnchorPairs on all edges are adjacent.
void TransitionGraph::check() const
{
    const TransitionGraph& transitionGraph = *this;

    BGL_FORALL_EDGES(e, transitionGraph, TransitionGraph) {

        const vertex_descriptor v0 = source(e, transitionGraph);
        const vertex_descriptor v1 = target(e, transitionGraph);

        const TransitionGraphVertex& vertex0 = transitionGraph[v0];
        const TransitionGraphVertex& vertex1 = transitionGraph[v1];

        SHASTA_ASSERT(vertex0.anchorPair.anchorIdB == vertex1.anchorPair.anchorIdA);
    }
}
