// Shasta2.
#include "AssemblyGraph.hpp"
#include "GTest.hpp"
#include "Options.hpp"
#include "performanceLog.hpp"
#include "TangleMatrix1.hpp"
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
        TangleMatrix1 tangleMatrix(assemblyGraph, inEdges, outEdges, html);
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
    performanceLog << timestamp << "AssemblyGraph::detangle begins." << endl;

    // Find the bubbles, including "haploid" bubbles.
    vector<Bubble> bubbles;
    findBubbles(bubbles, true);

    performanceLog << timestamp << "AssemblyGraph::detangle ends." << endl;
}
