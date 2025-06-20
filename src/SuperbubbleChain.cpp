#include "SuperbubbleChain.hpp"
#include "AssemblerOptions.hpp"
#include "deduplicate.hpp"
#include "PhasingGraph.hpp"
#include "Tangle.hpp"
#include "TangleMatrix.hpp"
using namespace shasta;

#include "algorithm.hpp"
#include "fstream.hpp"



void SuperbubbleChain::phase(
    AssemblyGraph& assemblyGraph,
    uint64_t superbubbleChainId)
{
    SuperbubbleChain& superbubbleChain = *this;
    const bool debug = false;

    const uint64_t m = 6;

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
    // non-trivial bubble. Each pair that can be phased
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
            if(debug) {
                cout << "Working on bubbles at positions " << position0 << " " << position1 << endl;
            }


            // Create a TangleMatrix between these two bubbles.
            TangleMatrix tangleMatrix(assemblyGraph,
                bubble0.internalEdges,
                bubble1.internalEdges,
                assemblyGraph.assemblerOptions.aDrift,
                assemblyGraph.assemblerOptions.bDrift);

            if(debug) {
                cout << "Tangle matrix:" << endl;
                for(uint64_t iEntrance=0; iEntrance<tangleMatrix.entrances.size(); iEntrance++) {
                    for(uint64_t iExit=0; iExit<tangleMatrix.exits.size(); iExit++) {
                        cout << tangleMatrix.tangleMatrix[iEntrance][iExit].orientedReadIds.size() << " ";
                    }
                    cout << endl;
                }
            }

            // Do the likelyHood ratio test (G test).
            const bool success = tangleMatrix.gTest(assemblyGraph.assemblerOptions.detangleEpsilon);

            if(not success) {
                if(debug) {
                    cout << "Likelihood ratio test was not successful." << endl;
                }
                continue;
            }

            if(debug) {
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

            if( tangleMatrix.hypotheses.front().isForwardInjective() or
                tangleMatrix.hypotheses.front().isBackwardInjective()) {
                phasingGraph.addEdge(position0, position1, tangleMatrix.hypotheses.front());

                if(debug) {
                    cout << "Added edge " << position0 << " " << position1 << endl;
                    cout << "Connectivity matrix for the best hypothesis:" << endl;
                    for(const auto& row: tangleMatrix.hypotheses.front().connectivityMatrix) {
                        for(const bool c: row) {
                            cout << int(c);
                        }
                        cout << endl;
                    }
                }
            }
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

    // Find the longest path for each connected component.
    phasingGraph.findLongestPaths();

    if(debug) {
        phasingGraph.writeGraphviz("PhasingGraph-" + to_string(superbubbleChainId) + ".dot");
    }



    // Loop over the longest paths. Each edge of a longest path generates a Tangle that can
    // be detangled using the Hypothesis stored in the edge and its connectivity matrix.
    std::set<AssemblyGraph::vertex_descriptor> removedVertices;
    for(uint64_t componentId=0; componentId<phasingGraph.longestPaths.size(); componentId++) {
        const vector<PhasingGraph::edge_descriptor>& longestPath = phasingGraph.longestPaths[componentId];
        for(const PhasingGraph::edge_descriptor e: longestPath) {
            const PhasingGraph::vertex_descriptor v0 = source(e, phasingGraph);
            const PhasingGraph::vertex_descriptor v1 = target(e, phasingGraph);

            const uint64_t position0 = phasingGraph[v0].position;
            const uint64_t position1 = phasingGraph[v1].position;

            const Superbubble& bubble0 = superbubbleChain[position0];
            const Superbubble& bubble1 = superbubbleChain[position1];

            if(debug) {
                cout << "Detangling between bubbles at positions " << position0 << " " << position1 << endl;
                cout << "Bubble at position " << position0 << ":";
                for(const AssemblyGraph::edge_descriptor e: bubble0.internalEdges) {
                    cout << " " << assemblyGraph[e].id;
                }
                cout << endl;
                cout << "Bubble at position " << position1 << ":";
                for(const AssemblyGraph::edge_descriptor e: bubble1.internalEdges) {
                    cout << " " << assemblyGraph[e].id;
                }
                cout << endl;
            }

            // Create a Tangle consisting of:
            // - The target vertex of bubble0.
            // - The source vertex of bubble1.
            // - All the  vertices of any bubbles or superbubbles in between.
            vector<AssemblyGraph::vertex_descriptor> tangleVertices;
            tangleVertices.push_back(bubble0.targetVertex);
            tangleVertices.push_back(bubble1.sourceVertex);
            for(uint64_t position=position0+1; position<position1; position++) {
                const Superbubble& superbubble = superbubbleChain[position];
                tangleVertices.push_back(superbubble.sourceVertex);
                tangleVertices.push_back(superbubble.targetVertex);
                copy(superbubble.internalVertices.begin(), superbubble.internalVertices.end(),
                    back_inserter(tangleVertices));
            }
            deduplicate(tangleVertices);

            // If any of these vertices have been removed, don't do anything.
            // This can happen occasionally.
            bool vertexWasRemoved = false;
            for(const AssemblyGraph::vertex_descriptor v: tangleVertices) {
                if(removedVertices.contains(v)) {
                    vertexWasRemoved = true;
                    break;
                }
            }
            if(vertexWasRemoved) {
                continue;
            }

            // Create the tangle from these vertices.
            Tangle tangle(assemblyGraph, tangleVertices,
                assemblyGraph.assemblerOptions.aDrift,
                assemblyGraph.assemblerOptions.bDrift);

            if(debug) {
                cout << "Tangle matrix:" << endl;
                for(uint64_t iEntrance=0; iEntrance<tangle.tangleMatrix->entrances.size(); iEntrance++) {
                    for(uint64_t iExit=0; iExit<tangle.tangleMatrix->exits.size(); iExit++) {
                        cout << tangle.tangleMatrix->tangleMatrix[iEntrance][iExit].orientedReadIds.size() << " ";
                    }
                    cout << endl;
                }
            }



            // Use the best hypothesis on this edge to detangle it.
            const TangleMatrix::Hypothesis& bestHypothesis = phasingGraph[e].bestHypothesis;
            const vector< vector<bool> >& connectivityMatrix = bestHypothesis.connectivityMatrix;
            if(debug) {
                cout << "Connectivity matrix for the best hypothesis:" << endl;
                for(uint64_t iEntrance=0; iEntrance<tangle.tangleMatrix->entrances.size(); iEntrance++) {
                    for(uint64_t iExit=0; iExit<tangle.tangleMatrix->exits.size(); iExit++) {
                        cout << int(connectivityMatrix[iEntrance][iExit]);
                    }
                    cout << endl;
                }
            }
            for(uint64_t iEntrance=0; iEntrance<tangle.tangleMatrix->entrances.size(); iEntrance++) {
                for(uint64_t iExit=0; iExit<tangle.tangleMatrix->exits.size(); iExit++) {
                    if(connectivityMatrix[iEntrance][iExit]) {
                        if(debug) {
                            cout << "Calling Tangle:: connect " << iEntrance << " " << iExit <<
                                ", anchor pair coverage " <<
                                tangle.tangleMatrix->tangleMatrix[iEntrance][iExit].orientedReadIds.size() << endl;
                        }
                        tangle.connect(iEntrance, iExit);
                    }
                }
            }
            if(debug) {
                cout << "Again, Tangle matrix:" << endl;
                for(uint64_t iEntrance=0; iEntrance<tangle.tangleMatrix->entrances.size(); iEntrance++) {
                    for(uint64_t iExit=0; iExit<tangle.tangleMatrix->exits.size(); iExit++) {
                        cout << tangle.tangleMatrix->tangleMatrix[iEntrance][iExit].orientedReadIds.size() << " ";
                    }
                    cout << endl;
                }
            }
            tangle.detangle();

            // Keep track of the vertices that were removed.
            for(vertex_descriptor v: tangle.removedVertices) {
                removedVertices.insert(v);
            }

        }
    }

}
