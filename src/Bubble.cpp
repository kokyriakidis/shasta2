#include "Bubble.hpp"
#include "AssemblyGraph.hpp"
#include "findConvergingVertex.hpp"
#include "Options.hpp"
#include "performanceLog.hpp"
#include "timestamp.hpp"
using namespace shasta2;



// Find Bubbles.
// The edges of each Bubble are sorted by id.
void AssemblyGraph::findBubbles(vector<Bubble>& bubbles, bool allowHaploid) const
{
    performanceLog << timestamp << "AssemblyGraph::findBubbles begins." << endl;

    const AssemblyGraph& assemblyGraph = *this;
    bubbles.clear();

    // Look at bubbles with source v0.
    std::map<vertex_descriptor, vector<edge_descriptor> > m;
    BGL_FORALL_VERTICES(v0, assemblyGraph, AssemblyGraph) {
        m.clear();
        BGL_FORALL_OUTEDGES(v0, e, assemblyGraph, AssemblyGraph) {
            const vertex_descriptor v1 = target(e, assemblyGraph);
            m[v1].push_back(e);
        }

        for(const auto& p: m) {
            const vertex_descriptor v1 = p.first;
            const vector<edge_descriptor>& edges = p.second;
            if(allowHaploid or (edges.size() > 1)) {
                Bubble bubble;
                bubble.v0 = v0;
                bubble.v1 = v1;
                bubble.edges = edges;
                std::ranges::sort(bubble.edges, orderById);
                bubbles.push_back(bubble);
            }
        }
    }

#if 0
    cout << "Found " << bubbles.size() << " bubbles." << endl;

    // Count the bubbles for each ploidy.
    vector<uint64_t> histogram;
    for(const Bubble& bubble: bubbles) {
        const uint64_t ploidy = bubble.edges.size();
        if(ploidy >= histogram.size()) {
            histogram.resize(ploidy+1, 0);
        }
        ++histogram[ploidy];
    }
    for(uint64_t ploidy=0; ploidy<histogram.size(); ploidy++) {
        const uint64_t frequency = histogram[ploidy];
        if(frequency) {
            cout << "Ploidy " << ploidy << ": " << frequency << " bubbles." << endl;
        }
    }
#endif

    performanceLog << timestamp << "AssemblyGraph::findBubbles ends." << endl;
}



// Find pairs of reverse complemented bubbles.
// The edges of the first bubble in each pair are sorted by id.
// The edges of the second bubble in each pair are sorted
// consistently with the ones in the first pair,
// that is, the reverse complement of p.first.edges[i] is p.second.edges[i].
void AssemblyGraph::findBubblePairs(vector<BubblePair>& bubblePairs) const
{
    performanceLog << timestamp << "AssemblyGraph::findBubblePairs begins." << endl;

    const AssemblyGraph& assemblyGraph = *this;
    bubblePairs.clear();

    // Each iteration of the loop generates a bubble and its reverse
    // complement. To generate each bubble only once we need to
    // keep track of the bubbles we already found.
    std::set< pair<vertex_descriptor, vertex_descriptor> > bubblesFound;

    // Look at bubbles with source v0.
    std::map<vertex_descriptor, vector<edge_descriptor> > m;
    BGL_FORALL_VERTICES(v0, assemblyGraph, AssemblyGraph) {

        // Gather edges, grouped by their target vertex.
        m.clear();
        BGL_FORALL_OUTEDGES(v0, e, assemblyGraph, AssemblyGraph) {
            const vertex_descriptor v1 = target(e, assemblyGraph);
            m[v1].push_back(e);
        }


        // Loop over groups with the same target vertex.
        // Each group with more than one edge generates
        // a pair of reverse complemented bubbles.
        for(auto& p: m) {
            const vertex_descriptor v1 = p.first;
            vector<edge_descriptor>& edges = p.second;
            if(edges.size() > 1) {

                // If we aleady generated this bubble, don't add it again.
                const vertex_descriptor v0Rc = assemblyGraph[v0].vRc;
                const vertex_descriptor v1Rc = assemblyGraph[v1].vRc;
                if(bubblesFound.contains(make_pair(v0, v1))) {
                    SHASTA2_ASSERT(bubblesFound.contains(make_pair(v1Rc, v0Rc)));
                    continue;
                }
                if(bubblesFound.contains(make_pair(v1Rc, v0Rc))) {
                    SHASTA2_ASSERT(bubblesFound.contains(make_pair(v0, v1)));
                    continue;
                }

                // Ok, if getting here we are going to generate a new pair of bubbles.
                bubblesFound.insert(make_pair(v0, v1));
                bubblesFound.insert(make_pair(v1Rc, v0Rc));

                std::ranges::sort(edges, orderById);
                auto& [bubble, bubbleRc] = bubblePairs.emplace_back();
                bubble.v0 = v0;
                bubble.v1 = v1;
                bubble.edges = edges;
                bubbleRc.v0 = v1Rc;
                bubbleRc.v1 = v0Rc;
                for(const edge_descriptor e: edges) {
                    bubbleRc.edges.push_back(assemblyGraph[e].eRc);
                }
            }
        }
    }

    performanceLog << timestamp << "AssemblyGraph::findBubblePairs ends." << endl;

}



// This finds all Superbubbles seen using the specified maxDistance.
// Some pairs of Superbubble can intersect (that is, they can have common edges).
void AssemblyGraph::findSuperbubbles(
    vector<Superbubble>& superbubbles) const
{
    const AssemblyGraph& assemblyGraph = *this;

    vector< pair<vertex_descriptor, vertex_descriptor> > forwardPairs;
    BGL_FORALL_VERTICES(vA, assemblyGraph, AssemblyGraph) {
        if(hasSelfEdge(vA)) {
            continue;
        }
        const vertex_descriptor vB = findConvergingVertexGeneral(assemblyGraph, vA, options.findSuperbubblesMaxDistance);
        if(vB != null_vertex()) {
            if(not hasSelfEdge(vB)) {
                forwardPairs.emplace_back(vA, vB);
                // cout << assemblyGraph[vA].id << "..." << assemblyGraph[vB].id << endl;
            }
        }
    }
    sort(forwardPairs.begin(), forwardPairs.end());
    // cout << "Found " << forwardPairs.size() << " forward pairs." << endl;

    const boost::reverse_graph<AssemblyGraph> reverseAssemblyGraph(assemblyGraph);
    vector< pair<vertex_descriptor, vertex_descriptor> > backwardPairs;
    BGL_FORALL_VERTICES(vA, assemblyGraph, AssemblyGraph) {
        if(hasSelfEdge(vA)) {
            continue;
        }
        const vertex_descriptor vB = findConvergingVertexGeneral(reverseAssemblyGraph, vA, options.findSuperbubblesMaxDistance);
        if(vB != null_vertex()) {
            if(not hasSelfEdge(vB)) {
                backwardPairs.emplace_back(vB, vA);
                // cout << assemblyGraph[vA].id << "..." << assemblyGraph[vB].id << endl;
            }
        }
    }
    sort(backwardPairs.begin(), backwardPairs.end());
    // cout << "Found " << backwardPairs.size() << " backward pairs." << endl;

    // Find pairs that appeared in both directions.
    vector< pair<vertex_descriptor, vertex_descriptor> > bidirectionalPairs;
    std::set_intersection(
        forwardPairs.begin(), forwardPairs.end(),
        backwardPairs.begin(), backwardPairs.end(),
        back_inserter(bidirectionalPairs));

    // Each of these pairs generates a Superbubble.
    superbubbles.clear();
    for(const auto& p: bidirectionalPairs) {
        superbubbles.emplace_back(assemblyGraph, p.first, p.second);

    }

}



void AssemblyGraph::findSuperbubbleChains(
    const vector<Superbubble>& superbubbles,
    vector<SuperbubbleChain>& superbubbleChains
    ) const
{
    // Index the superbubbles by their source and target vertex.
    std::map<vertex_descriptor, vector<uint64_t>, OrderById> mapBySource(orderById);
    std::map<vertex_descriptor, vector<uint64_t>, OrderById> mapByTarget(orderById);
    for(uint64_t superbubbleId=0; superbubbleId<superbubbles.size(); superbubbleId++) {
        const Superbubble& superbubble = superbubbles[superbubbleId];
        mapBySource[superbubble.sourceVertex].push_back(superbubbleId);
        mapByTarget[superbubble.targetVertex].push_back(superbubbleId);
    }

    // Sanity check: a vertex can only be a source or target of a single Superbubble.
    for(const auto& p: mapBySource) {
        SHASTA2_ASSERT(p.second.size() == 1);
    }
    for(const auto& p: mapByTarget) {
        SHASTA2_ASSERT(p.second.size() == 1);
    }

    // A vector to keep track of the Superbubbles we already added to a chain.
    vector<bool> wasUsed(superbubbles.size(), false);

    // Work vectors used below to construct chains.
    vector<uint64_t> forward;
    vector<uint64_t> backward;

    // Generate Superbubble chains.
    superbubbleChains.clear();
    for(uint64_t superbubbleId=0; superbubbleId<superbubbles.size(); superbubbleId++) {
        if(wasUsed[superbubbleId]) {
            continue;
        }
        // cout << "Starting a new chain at Superbubble " << superbubbleId << endl;

        // This Superbubble has not yet been added to any chain.
        // We will start a new chain here.

        // Create the forward portion of this chain.
        forward.clear();
        vertex_descriptor v = superbubbles[superbubbleId].targetVertex;
        while(true) {
            const auto it = mapBySource.find(v);
            if(it == mapBySource.end()) {
                break;
            }
            const vector<uint64_t>& nextVector = it->second;
            SHASTA2_ASSERT(nextVector.size() == 1);
            const uint64_t nextSuperbubbleId = nextVector.front();
            forward.push_back(nextSuperbubbleId);
            // cout << "Forward: " << nextSuperbubbleId << endl;

            v = superbubbles[nextSuperbubbleId].targetVertex;
        }

        // Create the backward portion of this chain.
        backward.clear();
        v = superbubbles[superbubbleId].sourceVertex;
        while(true) {
            const auto it = mapByTarget.find(v);
            if(it == mapByTarget.end()) {
                break;
            }
            const vector<uint64_t>& previousVector = it->second;
            SHASTA2_ASSERT(previousVector.size() == 1);
            const uint64_t previousSuperbubbleId = previousVector.front();
            backward.push_back(previousSuperbubbleId);
            // cout << "Backward: " << previousSuperbubbleId << endl;

            v = superbubbles[previousSuperbubbleId].sourceVertex;
        }

        // Now we can create the new SuperbubbleChain.
        superbubbleChains.emplace_back();
        SuperbubbleChain& superbubbleChain = superbubbleChains.back();
        std::reverse(backward.begin(), backward.end());
        for(const uint64_t id: backward) {
            wasUsed[id] = true;
            superbubbleChain.push_back(superbubbles[id]);
        }
        superbubbleChain.push_back(superbubbles[superbubbleId]);
        wasUsed[superbubbleId] = true;
        for(const uint64_t id: forward) {
            wasUsed[id] = true;
            superbubbleChain.push_back(superbubbles[id]);
        }
    }
}

