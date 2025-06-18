#include "SuperbubbleChain.hpp"
#include "AssemblerOptions.hpp"
#include "TangleMatrix.hpp"
using namespace shasta;

#include "algorithm.hpp"
#include "fstream.hpp"



void SuperbubbleChain::phase(
    AssemblyGraph& assemblyGraph,
    uint64_t superbubbleChainId)
{
    const bool debug = true;
    const uint64_t m = 3;

    if(debug) {
        cout << "Phasing superbubble chain " << superbubbleChainId << endl;
    }


    // An EdgeSet consists of one of the following:
    // - The internal edges of a non-trivial bubble.
    // - The source edges of a superbubble.
    // - The target edges of a superbubble.
    // To make sure all EdgeSets are disjoint,
    // we don't generate EdgeSets for superbubbles in which
    // some source edges are also target edges.
    // We will phase EdgeSets relative to each other.
    class EdgeSet {
    public:

        // Ths position in the SuperbubbleChain of the Superbubble (can be a Bubble)
        // that generated this edge set.
        uint64_t position;

        bool isInSuperbubble;
        bool isSuperbubbleTarget;

        vector<edge_descriptor> edges;
    };



    // Generate the EdgeSets.
    vector<EdgeSet> edgeSets;
    for(uint64_t position=0; position<size(); position++) {
        Superbubble& superbubble = (*this)[position];

        if(superbubble.isBubble()) {
            if(not superbubble.isTrivial()) {
                edgeSets.emplace_back();
                EdgeSet& edgeSet = edgeSets.back();
                edgeSet.position = position;
                edgeSet.isInSuperbubble = false;
                edgeSet.isSuperbubbleTarget = false;
                edgeSet.edges = superbubble.internalEdges;
            }

        } else {

            // Check if the source and target edges intersect.
            bool edgesIntersect = false;
            for(const edge_descriptor e: superbubble.sourceEdges) {
                if(target(e, assemblyGraph) == superbubble.targetVertex) {
                    edgesIntersect = true;
                    break;
                }
            }

            if(not edgesIntersect) {

                // Generate an EdgeSet for the source edges of this Superbubble.
                edgeSets.emplace_back();
                EdgeSet& sourceEdgeSet = edgeSets.back();
                sourceEdgeSet.position = position;
                sourceEdgeSet.isInSuperbubble = true;
                sourceEdgeSet.isSuperbubbleTarget = false;
                sourceEdgeSet.edges = superbubble.sourceEdges;

                // Generate an EdgeSet for the target edges of this Superbubble.
                edgeSets.emplace_back();
                EdgeSet& targetEdgeSet = edgeSets.back();
                targetEdgeSet.position = position;
                targetEdgeSet.isInSuperbubble = true;
                targetEdgeSet.isSuperbubbleTarget = true;
                targetEdgeSet.edges = superbubble.targetEdges;

            }
        }
    }

    // Write the EdgeSets.
    if(debug) {
        for(uint64_t edgeSetId=0; edgeSetId<edgeSets.size(); edgeSetId++) {
            const EdgeSet& edgeSet = edgeSets[edgeSetId];

            cout << "EdgeSet " << edgeSetId << ": ";
            if(edgeSet.isInSuperbubble) {
                if(edgeSet.isSuperbubbleTarget) {
                    cout << "Target edges of superbubble at position " << edgeSet.position << ":";
                } else {
                    cout << "Source edges of superbubble at position " << edgeSet.position << ":";
                }
            } else {
                cout << "Internal edges of bubble at position " << edgeSet.position << ":";
            }

            for(const edge_descriptor e: edgeSet.edges) {
                cout << " " << assemblyGraph[e].id;
            }
            cout << endl;
        }
    }





    // Phase pairs of EdgeSets.
    // Only consider pairs of edge sets whose index differs by no more than m.
    ofstream dot("PhasingGraph.dot");
    dot << "digraph PhasingGraph {" << endl;
    for(uint64_t index0=0; index0<edgeSets.size(); index0++) {
        const EdgeSet& edgeSet0 = edgeSets[index0];

        for(uint64_t index1=index0+1; index1<min(index0+m+1, edgeSets.size()); index1++) {
            const EdgeSet& edgeSet1 = edgeSets[index1];

            // Create a TangleMatrix between these two edge sets.
            TangleMatrix tangleMatrix(
                assemblyGraph,
                edgeSet0.edges,
                edgeSet1.edges,
                assemblyGraph.assemblerOptions.aDrift,
                assemblyGraph.assemblerOptions.bDrift);
            tangleMatrix.gTest(assemblyGraph.assemblerOptions.detangleEpsilon);

            SHASTA_ASSERT(not tangleMatrix.hypotheses.empty());
            const auto& bestHypothesis = tangleMatrix.hypotheses.front();

            // Figure out if this pair gives us usable information.
            const double bestG = bestHypothesis.G;
            double secondBestG = std::numeric_limits<double>::max();
            if(tangleMatrix.hypotheses.size() > 1) {
                secondBestG = tangleMatrix.hypotheses[1].G;
            }
            const double deltaLogP = secondBestG - bestG;
            bool isGood =
                (bestG <= assemblyGraph.assemblerOptions.detangleMaxLogP) and
                (deltaLogP >= assemblyGraph.assemblerOptions.detangleMinLogPDelta);

            if(debug) {
                cout << "Edge sets " << index0 << " " << index1 << endl;
                cout << "Best hypothesis:" << endl;
                for(const auto& row: bestHypothesis.connectivityMatrix) {
                    for(const auto value: row) {
                        cout << value << " ";
                    }
                    cout << endl;
                }
                cout << "logP " << bestHypothesis.G;
                if(tangleMatrix.hypotheses.size() > 1) {
                    cout << ", delta loGP " << tangleMatrix.hypotheses[1].G - bestHypothesis.G;
                }
                cout << (isGood ? " good" : " bad");
                cout << endl;
            }

            if(isGood) {
                dot << index0 << "->" << index1 << " [label=\"" <<
                    std::fixed << std::setprecision(1) << deltaLogP;
                for(const auto& row: bestHypothesis.connectivityMatrix) {
                    dot << "\\n";
                    for(const auto value: row) {
                        dot << value << " ";
                    }
                }
                dot << "\"];" << endl;
            }

        }

    }
    dot << "}" << endl;

}
