#include "SuperbubbleChain.hpp"
#include "AssemblerOptions.hpp"
#include "PhasingGraph.hpp"
#include "TangleMatrix.hpp"
using namespace shasta;

#include "algorithm.hpp"
#include "fstream.hpp"



void SuperbubbleChain::phase(
    AssemblyGraph& assemblyGraph,
    uint64_t superbubbleChainId)
{
    SuperbubbleChain& superbubbleChain = *this;
    const bool debug = true;

    const uint64_t m = 2;

    if(debug) {
        cout << "Phasing superbubble chain " << superbubbleChainId << endl;
    }

    // Create the PhasingGraph and generate its vertices.
    PhasingGraph phasingGraph;
    for(uint64_t position=0; position<size(); position++) {
        const Superbubble& bubble = superbubbleChain[position];
        if(not bubble.isBubble() or bubble.isTrivial()) {
            continue;
        }
        phasingGraph.addVertex(position);
    }


    // Loop over close pairs of Superbubbles that consist of a single
    // non-trivial bubble. Eaqch pair that can be phased
    // generates an edge of the PhasingGraph.
    for(uint64_t position0=0; position0<size(); position0++) {
        const Superbubble& bubble0 = superbubbleChain[position0];
        if(not bubble0.isBubble() or bubble0.isTrivial()) {
            continue;
        }
        uint64_t n0 = 0;
        for(uint64_t position1=position0+1; position1<size(); position1++) {
            const Superbubble& bubble1 = superbubbleChain[position1];
            if(not bubble1.isBubble() or bubble1.isTrivial()) {
                continue;
            }
            if(n0 > m) {
                continue;
            }
            ++n0;
            if(false) {
                cout << "Working on bubbles at positions " << position0 << " " << position1 << endl;
            }


            // Create a TangleMatrix between these two bubbles.
            TangleMatrix tangleMatrix(assemblyGraph,
                bubble0.internalEdges,
                bubble1.internalEdges,
                assemblyGraph.assemblerOptions.aDrift,
                assemblyGraph.assemblerOptions.bDrift);

            // Do the Likely hood rato test (G test).
            const bool success = tangleMatrix.gTest(assemblyGraph.assemblerOptions.detangleEpsilon);

            if(not success) {
                if(debug) {
                    cout << "Likelihood ratio test was not successful." << endl;
                }
                continue;
            }

            if(false) {
                cout << "Best hypothesis " << tangleMatrix.hypotheses.front().G;
                if(tangleMatrix.hypotheses.size() > 1) {
                    cout << ", second best hypothesis " << tangleMatrix.hypotheses[1].G;
                }
                cout << endl;
            }

            // Check if the best hypothesis satisfies our options.
            const double bestG = tangleMatrix.hypotheses.front().G;
            if(bestG > assemblyGraph.assemblerOptions.detangleMaxLogP) {
                if(false) {
                    cout << "Best hypothesis G is too high." << endl;
                }
                continue;
            }
            if(tangleMatrix.hypotheses.size() > 1) {
                const double secondBestG = tangleMatrix.hypotheses[1].G;
                if(secondBestG - bestG < assemblyGraph.assemblerOptions.detangleMinLogPDelta) {
                    if(false) {
                        cout << "Second best hypothesis G is too low." << endl;
                    }
                    continue;
                }
            }

            phasingGraph.addEdge(position0, position1);
        }
    }

    if(debug) {
        cout << "This superbubble chain has " << num_vertices(phasingGraph) <<
            " non-trivial bubbles." << endl;
        cout << "The phasing graph has " << num_vertices(phasingGraph) <<
            " vertices and " << num_edges(phasingGraph) << " edges." << endl;
    }

    if(num_vertices(phasingGraph) < 2) {
        if(debug) {
            cout << "There is nothing to phase." << endl;
        }
        return;
    }

    // Remove isolated vertices.
    const uint64_t removedIsolatedVertexCount = phasingGraph.removeIsolatedVertices();
    if(debug) {
        cout << removedIsolatedVertexCount << " isolated vertices were removed." << endl;
    }
    if(num_vertices(phasingGraph) < 2) {
        return;
    }

    // Compute connected components.
    phasingGraph.computeConnectedComponents();
    if(debug) {
        if(phasingGraph.components.size() == 1) {
            cout << "The phasing graph has a single connected component of size " <<
                phasingGraph.components.front().size() << " vertices." << endl;
        } else {
            cout << "The phasing graph has " << phasingGraph.components.size() <<
                " connected components of sizes";
            for(const vector<uint64_t>& component: phasingGraph.components) {
                cout << " " << component.size();
            }
            cout << endl;
        }
    }

    if(debug) {
        phasingGraph.writeGraphviz("PhasingGraph-" + to_string(superbubbleChainId) + ".dot");
    }

}
