// Shasta2.
#include "AssemblyGraph.hpp"
#include "color.hpp"
#include "DisjointSets.hpp"
#include "GTest.hpp"
#include "Options.hpp"
#include "performanceLog.hpp"
#include "TangleGraph.hpp"
#include "TangleMatrix.hpp"
#include "timestamp.hpp"
using namespace shasta2;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>



// This only handles simple vertex tangles with a 2 by 2 tangle matrix.
void AssemblyGraph::detangleVertices()
{
    AssemblyGraph& assemblyGraph = *this;
    const bool debug = false;
    ostream html(0);

    const vector< vector<bool> > inPhaseConnectivityMatrix    = { {true , false}, {false, true } };
    const vector< vector<bool> > outOfPhaseConnectivityMatrix = { {false, true }, {true , false} };

    // For now we only handle vertices with in-degree and out-degree equal to 2.
    // Gather the ones with id less that the id of their reverse complement.
    vector<vertex_descriptor> candidateVertices;
    BGL_FORALL_VERTICES(v, assemblyGraph, AssemblyGraph) {
        if(in_degree(v, assemblyGraph) != 2) {
            continue;
        }
        if(out_degree(v, assemblyGraph) != 2) {
            continue;
        }
        const AssemblyGraphVertex& vertex = assemblyGraph[v];
        const vertex_descriptor vRc = vertex.vRc;
        SHASTA2_ASSERT(vRc != null_vertex());
        SHASTA2_ASSERT(vRc != v);
        const AssemblyGraphVertex& vertexRc = assemblyGraph[vRc];
        if(vertex.id < vertexRc.id) {
            candidateVertices.push_back(v);
        }
    }
    cout << "Found " << candidateVertices.size() << " candidate vertex pairs for detangling." << endl;



    // Now loop over out detangling candidates.
    vector<edge_descriptor> inEdges;
    vector<edge_descriptor> outEdges;
    vector<edge_descriptor> inEdgesRc;
    vector<edge_descriptor> outEdgesRc;
    for(const vertex_descriptor v: candidateVertices) {
        const AssemblyGraphVertex& vertex = assemblyGraph[v];
        const vertex_descriptor vRc = vertex.vRc;
        SHASTA2_ASSERT(vRc != null_vertex());
        SHASTA2_ASSERT(vRc != v);
        const AssemblyGraphVertex& vertexRc = assemblyGraph[vRc];

        if(debug) {
            cout << "Working on vertex pair " << vertex.id << " " << vertexRc.id << endl;
        }

        // Gather incoming/outgoing edges and sp we can run some sanity checks on them.
        inEdges.clear();
        outEdges.clear();
        inEdgesRc.clear();
        outEdgesRc.clear();
        BGL_FORALL_INEDGES(v, e, assemblyGraph, AssemblyGraph) {
            inEdges.push_back(e);
        }
        BGL_FORALL_OUTEDGES(v, e, assemblyGraph, AssemblyGraph) {
            outEdges.push_back(e);
        }
        BGL_FORALL_INEDGES(vRc, eRc, assemblyGraph, AssemblyGraph) {
            inEdgesRc.push_back(eRc);
        }
        BGL_FORALL_OUTEDGES(vRc, eRc, assemblyGraph, AssemblyGraph) {
            outEdgesRc.push_back(eRc);
        }
        std::ranges::sort(inEdges, orderById);
        std::ranges::sort(outEdges, orderById);
        std::ranges::sort(inEdgesRc, orderById);
        std::ranges::sort(outEdgesRc, orderById);

        if(debug) {
            cout << "inEdges:";
            for(const edge_descriptor e: inEdges) {
                cout << " " << assemblyGraph[e].id;
            }
            cout << endl;
            cout << "outEdges:";
            for(const edge_descriptor e: outEdges) {
                cout << " " << assemblyGraph[e].id;
            }
            cout << endl;
            cout << "inEdgesRc:";
            for(const edge_descriptor eRc: inEdgesRc) {
                cout << " " << assemblyGraph[eRc].id;
            }
            cout << endl;
            cout << "outEdgesRc:";
            for(const edge_descriptor eRc: outEdgesRc) {
                cout << " " << assemblyGraph[eRc].id;
            }
            cout << endl;
        }

        // If an edge is both an in-edge and out-edge, skip this vertex.
        bool hasCycle = false;
        for(const edge_descriptor e: inEdges) {
            if(std::ranges::contains(outEdges, e)) {
                hasCycle = true;
            }
        }
        for(const edge_descriptor e: outEdges) {
            if(std::ranges::contains(inEdges, e)) {
                hasCycle = true;
            }
        }
        for(const edge_descriptor eRc: inEdgesRc) {
            if(std::ranges::contains(outEdgesRc, eRc)) {
                hasCycle = true;
            }
        }
        for(const edge_descriptor eRc: outEdgesRc) {
            if(std::ranges::contains(inEdgesRc, eRc)) {
                hasCycle = true;
            }
        }
        if(hasCycle) {
            if(debug) {
                cout << "Skipped due to cycle." << endl;
            }
            continue;
        }



        // Check the reverse complements of the in-edges and out-edges.
        for(const edge_descriptor e: inEdges) {
            const edge_descriptor eRc = assemblyGraph[e].eRc;
            SHASTA2_ASSERT(std::ranges::contains(outEdgesRc, eRc));
        }
        for(const edge_descriptor e: outEdges) {
            const edge_descriptor eRc = assemblyGraph[e].eRc;
            SHASTA2_ASSERT(std::ranges::contains(inEdgesRc, eRc));
        }
        for(const edge_descriptor eRc: inEdgesRc) {
            const edge_descriptor e = assemblyGraph[eRc].eRc;
            SHASTA2_ASSERT(std::ranges::contains(outEdges, e));
        }
        for(const edge_descriptor eRc: outEdgesRc) {
            const edge_descriptor e = assemblyGraph[eRc].eRc;
            SHASTA2_ASSERT(std::ranges::contains(inEdges, e));
        }

        // Check that the two reverse complemented tangles are entirely disjoint
        // (they don't share any edges).
        for(const edge_descriptor e: inEdges) {
            SHASTA2_ASSERT(not std::ranges::contains(inEdgesRc, e));
            SHASTA2_ASSERT(not std::ranges::contains(outEdgesRc, e));
        }
        for(const edge_descriptor e: outEdges) {
            SHASTA2_ASSERT(not std::ranges::contains(inEdgesRc, e));
            SHASTA2_ASSERT(not std::ranges::contains(outEdgesRc, e));
        }

        // Create the tangle matrix.
        TangleMatrix tangleMatrix(assemblyGraph, inEdges, outEdges, html);
        if(debug) {
            cout << "Tangle  matrix: ";
            for(const auto& row: tangleMatrix.tangleMatrix) {
                std::ranges::copy(row, ostream_iterator<uint64_t>(cout, " "));
            }
            cout << endl;
        }

        // Run the G-test on this tangle matrix.
        GTest gTest(tangleMatrix.tangleMatrix, options.detangleEpsilon, false, false);

        // If the G-test failed, don't detangle.
        if(not gTest.success) {
            if(debug) {
                cout << "Likelihood ratio test was not successful." << endl;
            }
            continue;
        }

        const auto& bestHypothesis = gTest.hypotheses.front();
        const double bestG = bestHypothesis.G;
        if(debug) {
            cout << "Best hypothesis:";
            for(const auto& row: bestHypothesis.connectivityMatrix) {
                for(const bool b: row) {
                    cout << (b ? " 1" : " 0");
                }
            }
            cout << endl;
            cout << "G = " << bestG;
            if(gTest.hypotheses.size() > 1) {
                const double secondBestG = gTest.hypotheses[1].G;
                cout << ", second best G = " << secondBestG;
            }
            cout << endl;
        }

        // Check if the best hypothesis satisfies our options.
        if(bestG > options.detangleMaxLogP) {
            if(debug) {
                cout << "Best hypothesis G is too high." << endl;
            }
            continue;
        }
        if(gTest.hypotheses.size() > 1) {
            const double secondBestG = gTest.hypotheses[1].G;
            if(secondBestG - bestG < assemblyGraph.options.detangleMinLogPDelta) {
                if(debug) {
                    cout << "Second best hypothesis G is too low." << endl;
                }
                continue;
            }
        }

        // Figure out if in phase or out of phase.
        bool isInPhase = (bestHypothesis.connectivityMatrix == inPhaseConnectivityMatrix);
        bool isOutOfPhase = (bestHypothesis.connectivityMatrix == outOfPhaseConnectivityMatrix);
        if(not (isInPhase or isOutOfPhase)) {
            if(debug) {
                cout << "Not detangling because the best hypothesis is not a permutation." << endl;
            }
            continue;
        }
        if(debug) {
            if(isInPhase) {
                cout << "This vertex and its reverse complement will be detangled in phase." << endl;
            } else {
                SHASTA2_ASSERT(isOutOfPhase);
                cout << "This vertex and its reverse complement will be detangled out of phase." << endl;
            }
        }



        // Loop over the entrance/exit pairs to be connected.
        for(uint64_t i=0; i<2; i++) {
            const edge_descriptor entrance = inEdges[i];
            const edge_descriptor exit = (isInPhase ? outEdges[i] : outEdges[1-i]);

            // Combine this entrance and exit into a single edge.
            const vertex_descriptor v0 = source(entrance, assemblyGraph);
            const vertex_descriptor v1 = target(exit, assemblyGraph);
            auto[eNew, wasAdded] = add_edge(v0, v1, AssemblyGraphEdge(nextEdgeId++), assemblyGraph);
            SHASTA2_ASSERT(wasAdded);
            AssemblyGraphEdge& edgeNew = assemblyGraph[eNew];

            // Add the entrance steps to the new edge.
            edgeNew.swapSteps(assemblyGraph[entrance]);

            // Add the exit steps to the new edge.
            const vector<AssemblyGraphEdgeStep>& exitSteps = assemblyGraph[exit];
            std::ranges::copy(exitSteps, back_inserter(edgeNew));

            // Create the reverse complement of the new edge.
            createReverseComplementEdge(eNew);

            // Now we can remove the entrance, the exit, and their reverse complements.
            boost::remove_edge(assemblyGraph[entrance].eRc, assemblyGraph);
            boost::remove_edge(assemblyGraph[exit].eRc, assemblyGraph);
            boost::remove_edge(entrance, assemblyGraph);
            boost::remove_edge(exit, assemblyGraph);
        }

        // Now we can remove the vertex we detangled and its reverse complement.
        SHASTA2_ASSERT(in_degree(v, assemblyGraph) == 0);
        SHASTA2_ASSERT(out_degree(v, assemblyGraph) == 0);
        SHASTA2_ASSERT(in_degree(vRc, assemblyGraph) == 0);
        SHASTA2_ASSERT(out_degree(vRc, assemblyGraph) == 0);
        boost::remove_vertex(v, assemblyGraph);
        boost::remove_vertex(vRc, assemblyGraph);

    }

}



void AssemblyGraph::detangle()
{
    // EXPOSE WHEN CODE STABILIZES.
    const uint64_t lengthThreshold = 30000;

    performanceLog << timestamp << "AssemblyGraph::detangle begins." << endl;
    AssemblyGraph& assemblyGraph = *this;

    // Find the BubbleChains.
    vector<BubbleChain> bubbleChains;
    findBubbleChains(bubbleChains);
    writeBubbleChains("BubbleChains.csv", bubbleChains);
    writeBubbleChainsForBandage("BubbleChains-Bandage.csv", bubbleChains);

    // Store separately short and long BubbleChains.
    vector<BubbleChain> shortBubbleChains;
    vector<BubbleChain> longBubbleChains;
    for(const BubbleChain& bubbleChain: bubbleChains) {
        if(bubbleChain.maxLength(assemblyGraph) >= lengthThreshold) {
            longBubbleChains.emplace_back(bubbleChain);
        } else {
            shortBubbleChains.emplace_back(bubbleChain);
        }
    }
    writeBubbleChains("LongBubbleChains.csv", longBubbleChains);
    writeBubbleChainsForBandage("LongBubbleChains-Bandage.csv", longBubbleChains);
    cout << "Found " << longBubbleChains.size() << " long bubble chains." << endl;

    // Map the vertices to integers.
    // This is needed below to compute connected components.
    std::map<vertex_descriptor, uint64_t> vertexIndexMap;
    std::vector<vertex_descriptor> vertexTable;
    BGL_FORALL_VERTICES(v, assemblyGraph, AssemblyGraph) {
        vertexIndexMap.insert(make_pair(v, vertexTable.size()));
        vertexTable.push_back(v);
    }



    // Compute connected components using only edges that belong to short BubbleChains.
    // Each non-trivial connected component will become a tangle.
    DisjointSets disjointSets(vertexTable.size());
    for(const BubbleChain& bubbleChain: shortBubbleChains) {

        for(const Bubble& bubble: bubbleChain) {
            disjointSets.unionSet(vertexIndexMap.at(bubble.v0), vertexIndexMap.at(bubble.v1));
        }
    }

    // Get the connected components with two or more vertices.
    vector< vector<uint64_t> > componentsVertexIndexes;
    disjointSets.gatherComponents(2, componentsVertexIndexes);

    // Convert vertex indexes to vertex descriptors.
    // Sort each component so we can do binary searches.
    vector< vector<vertex_descriptor> > tangles;
    for(const auto& componentVertexIndexes: componentsVertexIndexes) {
        vector<vertex_descriptor>& tangle = tangles.emplace_back();
        for(const uint64_t vertexIndex: componentVertexIndexes) {
            tangle.push_back(vertexTable[vertexIndex]);
        }
        sort(tangle.begin(), tangle.end(), orderById);
    }


    // Write a csv file that can be imported into Bandage to see
    // the tangles.
    {
        ofstream csv("Tangles-Bandage.csv");
        csv << "Segment,Tangle,Color\n";
        for(uint64_t tangleId=0; tangleId<tangles.size(); tangleId++) {
            const string color = randomHslColor(tangleId, 0.75, 0.5);
            const vector<vertex_descriptor>& tangle = tangles[tangleId];
            for(const vertex_descriptor v0: tangle) {
                BGL_FORALL_OUTEDGES(v0, e, assemblyGraph, AssemblyGraph) {
                    const vertex_descriptor v1 = target(e, assemblyGraph);
                    if(binary_search(tangle.begin(), tangle.end(), v1, orderById)) {
                        csv << assemblyGraph[e].id << ",";
                        csv << tangleId << ",";
                        csv << color << "\n";
                    }
                }
            }
        }
    }


    // Create the TangleGraph and use it to detangle.
    TangleGraph tangleGraph(assemblyGraph, longBubbleChains, tangles);
    tangleGraph.detangle();

    performanceLog << timestamp << "AssemblyGraph::detangle ends." << endl;
}
