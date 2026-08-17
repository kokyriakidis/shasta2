// Shasta2.
#include "AssemblyGraph.hpp"
#include "color.hpp"
#include "deduplicate.hpp"
#include "DisjointSets.hpp"
#include "GTest.hpp"
#include "Options.hpp"
#include "performanceLog.hpp"
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
    for(uint64_t iteration=0; ; ++iteration) {
        cout << "Detangle iteration " << iteration << " begins." << endl;
        write("Before-Iteration-" + to_string(iteration));

        const bool somethingWasDone = detangleIteration();
        if(not somethingWasDone) {
            break;
        }
        write("After-Iteration-" + to_string(iteration));

        strandSymmetricCompress();
    }

    check();
}



bool AssemblyGraph::detangleIteration()
{
    // EXPOSE WHEN CODE STABILIZES.
    const uint64_t lengthThreshold = 30000;

    // cout << timestamp << "AssemblyGraph::detangleIteration begins." << endl;
    AssemblyGraph& assemblyGraph = *this;

    // Find the BubbleChains.
    vector<BubbleChain> bubbleChains;
    findBubbleChains(bubbleChains);
    writeBubbleChains("BubbleChains.csv", bubbleChains);
    writeBubbleChainsForBandage("BubbleChains-Bandage.csv", bubbleChains);

    // Ony keep the short BubbleChains.
    // These will be used to define our tangles.
    vector<BubbleChain> shortBubbleChains;
    for(const BubbleChain& bubbleChain: bubbleChains) {
        if(bubbleChain.maxLength(assemblyGraph) < lengthThreshold) {
            shortBubbleChains.emplace_back(bubbleChain);
        }
    }
    writeBubbleChains("ShortBubbleChains.csv", shortBubbleChains);
    writeBubbleChainsForBandage("ShortBubbleChains-Bandage.csv", shortBubbleChains);
    cout << "Found " << shortBubbleChains.size() << " short bubble chains." << endl;

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
    // These will be our tangles.
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
        std::ranges::sort(tangle, orderById);
    }
    cout << "Found " << tangles.size() << " tangles." << endl;

    // Create a map that gives the tangle each vertex belongs to, if any.
    std::map<vertex_descriptor, uint64_t> tangleMap;
    for(uint64_t tangleId=0; tangleId<tangles.size(); tangleId++) {
        const vector<vertex_descriptor>& tangle = tangles[tangleId];
        for(const vertex_descriptor v: tangle) {
            tangleMap.insert(make_pair(v, tangleId));
        }
    }

    // Find the reverse complement of each tangle.
    vector<uint64_t> tangleRc(tangles.size(), invalid<uint64_t>);
    for(uint64_t tangleId=0; tangleId<tangles.size(); tangleId++) {
        const vector<vertex_descriptor>& tangle = tangles[tangleId];
        const vertex_descriptor v = tangle.front();
        const vertex_descriptor vRc = assemblyGraph[v].vRc;
        tangleRc[tangleId] = tangleMap.at(vRc);
    }

    // Write a csv file that can be imported into Bandage to see
    // the tangles.
    {
        ofstream csv("Tangles-Bandage.csv");
        csv << "Segment,Tangle,TangleRc,Color\n";
        for(uint64_t tangleId=0; tangleId<tangles.size(); tangleId++) {
            const string color = randomHslColor(tangleId, 0.75, 0.5);
            const vector<vertex_descriptor>& tangle = tangles[tangleId];
            for(const vertex_descriptor v0: tangle) {
                BGL_FORALL_OUTEDGES(v0, e, assemblyGraph, AssemblyGraph) {
                    const vertex_descriptor v1 = target(e, assemblyGraph);
                    if(binary_search(tangle.begin(), tangle.end(), v1, orderById)) {
                        csv << id(e) << ",";
                        csv << tangleId << ",";
                        csv << tangleRc[tangleId] << ",";
                        csv << color << "\n";
                    }
                }
            }
        }
    }

    // Detangle each pair of reverse complemented tangles.
    bool somethingWasDone = false;
    for(uint64_t tangleId=0; tangleId<tangles.size(); tangleId++) {
        if(tangleId <= tangleRc[tangleId]) {
            const vector<vertex_descriptor>& tangle = tangles[tangleId];
            const bool success = detangleStrandSymmetric(tangleId, tangle);
            somethingWasDone = somethingWasDone or success;
        }
    }

    /*
    cout << timestamp << "AssemblyGraph::detangleIteration ends. " <<
        (somethingWasDone ? "Something" : "Nothing") << " was done." << endl;
    */

    return somethingWasDone;
}



bool AssemblyGraph::detangleStrandSymmetric(
    uint64_t tangleId,
    const vector<vertex_descriptor>& tangleVertices)
{
    AssemblyGraph& assemblyGraph = *this;

    const bool debug = false;
    if(debug) {
        cout << "Working on tangle " << tangleId <<
            " with " << tangleVertices.size() << " vertices." <<endl;
    }
    SHASTA2_ASSERT(std::ranges::is_sorted(tangleVertices, orderById));

    // Figure out if this tangle is self-complementary.
    const vertex_descriptor v = tangleVertices.front();
    const vertex_descriptor vRc = assemblyGraph[v].vRc;
    const bool isSelfComplementary = (std::ranges::binary_search(tangleVertices, vRc, orderById));
    if(debug) {
        cout << "This tangle is" << (isSelfComplementary ? "" : " not") <<
            " self-complementary." << endl;
    }

    // For now we don't handle the self-complementary case.
    if(isSelfComplementary) {
        if(debug) {
            cout << "Detangling for self-complementary tangles is not implemented." << endl;
        }
        return false;
    }

    // Get the entrances.
    vector<Segment> entrances;
    for(const vertex_descriptor v0: tangleVertices) {
        BGL_FORALL_INEDGES(v0, e, assemblyGraph, AssemblyGraph) {
            const vertex_descriptor v1 = source(e, assemblyGraph);
            if(not std::ranges::binary_search(tangleVertices, v1, orderById)) {
                entrances.push_back(e);
            }
        }
    }
    std::ranges::sort(entrances, orderById);

    // Get the exits.
    vector<Segment> exits;
    for(const vertex_descriptor v0: tangleVertices) {
        BGL_FORALL_OUTEDGES(v0, e, assemblyGraph, AssemblyGraph) {
            const vertex_descriptor v1 = target(e, assemblyGraph);
            if(not std::ranges::binary_search(tangleVertices, v1, orderById)) {
                exits.push_back(e);
            }
        }
    }
    std::ranges::sort(exits, orderById);

    if(debug) {
        cout << "This tangle has " << entrances.size() << " entrances:";
        for(const Segment entrance: entrances) {
            cout << " " << id(entrance);
        }
        cout << endl;
        cout << "This tangle has " << exits.size() << " exits:";
        for(const Segment exit: exits) {
            cout << " " << id(exit);
        }
        cout << endl;
    }

    if(entrances.empty() or exits.empty()) {
        if(debug) {
            cout << "Not detangling because there are no entrances or no exits." << endl;
        }
        return false;
    }

    // If any entrances are also exit, don't detangle for now.
    bool hasEntranceExit = false;
    for(const Segment entrance: entrances) {
        if(std::ranges::contains(exits, entrance)) {
            hasEntranceExit = true;
        }
    }
    if(hasEntranceExit) {
        if(debug) {
            cout << "Not detangling because an entrance is also an exit."<< endl;
        }
        return false;
    }


    // Create the TangleMatrix.
    ostream noOutput(0);
    const TangleMatrix tangleMatrix(assemblyGraph, entrances, exits, noOutput);
    if(debug) {
        cout << "Tangle matrix:" << endl;
        cout << ",";
        for(const Segment exit: exits) {
            cout << id(exit) << ",";
        }
        cout << endl;
        for(uint64_t i=0; i<entrances.size(); i++) {
            cout << id(entrances[i]) << ",";
            for(uint64_t j=0; j<exits.size(); j++) {
                cout << tangleMatrix.tangleMatrix[i][j] << ",";
            }
            cout << endl;
        }
    }

    // Run the G-test on this tangle matrix.
    GTest gTest(tangleMatrix.tangleMatrix, assemblyGraph.options.detangleEpsilon, false, false);

    // If the G-test failed, don't detangle.
    if(not gTest.success) {
        if(debug) {
            cout << "Likelihood ratio test was not successful." << endl;
        }
        return false;
    }

    const auto& bestHypothesis = gTest.hypotheses.front();
    const double bestG = bestHypothesis.G;
    if(debug) {
        cout << "Best hypothesis:" << endl;
        cout << ",";
        for(const AssemblyGraph::edge_descriptor exit: exits) {
            cout << assemblyGraph[exit].id << ",";
        }
        cout << endl;
        for(uint64_t i=0; i<entrances.size(); i++) {
            cout << assemblyGraph[entrances[i]].id << ",";
            for(uint64_t j=0; j<exits.size(); j++) {
                cout << int(bestHypothesis.connectivityMatrix[i][j]) << ",";
            }
            cout << endl;
        }
        cout << "G = " << bestG;
        if(gTest.hypotheses.size() > 1) {
            const double secondBestG = gTest.hypotheses[1].G;
            cout << ", second best G = " << secondBestG;
        }
        cout << endl;
    }

    // Check if the best hypothesis satisfies our options.
    if(bestG > assemblyGraph.options.detangleMaxLogP) {
        if(debug) {
            cout << "Best hypothesis G is too high." << endl;
        }
        return false;
    }
    if(gTest.hypotheses.size() > 1) {
        const double secondBestG = gTest.hypotheses[1].G;
        if(secondBestG - bestG < assemblyGraph.options.detangleMinLogPDelta) {
            if(debug) {
                cout << "Second best hypothesis G is too low." << endl;
            }
            return false;
        }
    }

    // Check if  we can connect the entrance/exit pairs
    // described by the connectivity matrix for the best hypothesis.
    for(uint64_t i=0; i<entrances.size(); i++) {
        const Segment entrance = entrances[i];
        for(uint64_t j=0; j<exits.size(); j++) {
            if(bestHypothesis.connectivityMatrix[i][j]) {
                const Segment exit = exits[j];
                if(not assemblyGraph.canConnect(entrance, exit)) {
                    if(debug) {
                        cout << "Not detangling because can't connect " << id(entrance) <<
                            " with " << id(exit) << endl;
                    }
                    return false;
                }
            }
        }
    }

    if(debug) {
        if(isSelfComplementary) {
            cout << "This tangle will be detangled." << endl;
        } else {
            cout << "This tangle and its reverse complement will be detangled." << endl;
        }
    }

    // The code below relies on these two assumptions, checked above.
    SHASTA2_ASSERT(not isSelfComplementary);
    SHASTA2_ASSERT(not hasEntranceExit);



    // Remove all the Segments internal to this tangle.
    vector<Segment> edgesToBeRemoved;
    for(const vertex_descriptor v0: tangleVertices) {
        BGL_FORALL_OUTEDGES(v0, e, assemblyGraph, AssemblyGraph) {
            const AssemblyGraph::vertex_descriptor v1 = target(e, assemblyGraph);
            if(std::ranges::binary_search(tangleVertices, v1)) {
                edgesToBeRemoved.push_back(e);
                edgesToBeRemoved.push_back(assemblyGraph[e].eRc);
            }
        }
    }
    deduplicate(edgesToBeRemoved);
    for(const Segment e: edgesToBeRemoved) {
        boost::remove_edge(e, assemblyGraph);
    }



    // Create copies of the entrances disconnected at the end
    // and of the exits disconnected at the beginning.
    vector<Segment> newEntrances;
    for(const Segment entranceOld: entrances) {
        newEntrances.push_back(disconnectAtEnd(entranceOld));
    }
    vector<Segment> newExits;
    for(const Segment exitOld: exits) {
        newExits.push_back(disconnectAtBeginning(exitOld));
    }

    // Make the connections described by the connectivity matrix
    // of the best hypothesis.
    for(uint64_t i=0; i<newEntrances.size(); i++) {
        const AssemblyGraph::edge_descriptor newEntrance = newEntrances[i];
        for(uint64_t j=0; j<newExits.size(); j++) {
            if(bestHypothesis.connectivityMatrix[i][j]) {
                const AssemblyGraph::edge_descriptor newExit = newExits[j];
                const edge_descriptor eNew = connect(newEntrance, newExit);
                createReverseComplementEdge(eNew);
            }
        }
    }

    return true;
}

