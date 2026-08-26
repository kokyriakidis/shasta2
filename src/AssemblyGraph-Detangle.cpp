// Shasta2.
#include "AssemblyGraph.hpp"
#include "color.hpp"
#include "deduplicate.hpp"
#include "DisjointSets.hpp"
#include "GTest.hpp"
#include "html.hpp"
#include "Options.hpp"
#include "performanceLog.hpp"
#include "Tangle.hpp"
#include "TangleMatrix.hpp"
#include "timestamp.hpp"
using namespace shasta2;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>



void AssemblyGraph::detangleVertices()
{
    performanceLog << timestamp << "AssemblyGraph::detangleVertices begins." << endl;

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

        if(id(v) < id(vRc)) {
            candidateVertices.push_back(v);
        }
    }
    cout << "Found " << candidateVertices.size() << " candidate vertex pairs for detangling." << endl;



    // Now loop over out detangling candidates.
    vector<edge_descriptor> inEdges;        // in-edges of v
    vector<edge_descriptor> outEdges;       // out-edges of v
    vector<edge_descriptor> inEdgesRc;      // in-edges of vRc
    vector<edge_descriptor> outEdgesRc;     // out-edges of vRc
    for(const vertex_descriptor v: candidateVertices) {
        const AssemblyGraphVertex& vertex = assemblyGraph[v];
        const vertex_descriptor vRc = vertex.vRc;
        SHASTA2_ASSERT(vRc != null_vertex());
        SHASTA2_ASSERT(vRc != v);

        if(debug) {
            cout << "Working on vertex pair " << id(v) << " " << id(vRc) << endl;
        }

        // Gather incoming/outgoing edges.
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
                cout << " " << id(e);
            }
            cout << endl;
            cout << "outEdges:";
            for(const edge_descriptor e: outEdges) {
                cout << " " << id(e);
            }
            cout << endl;
            cout << "inEdgesRc:";
            for(const edge_descriptor e: inEdgesRc) {
                cout << " " << id(e);
            }
            cout << endl;
            cout << "outEdgesRc:";
            for(const edge_descriptor e: outEdgesRc) {
                cout << " " << id(e);
            }
            cout << endl;
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



        // If getting here, this vertex and its reverse complement will be detangled.
        // To do this, we need to disconnect (at the beginning or end) the Segments involved.
        // To do this we use disconnectAtBeginning and disconnectAtEnd, which change
        // edge_descriptors but leave the ids unchanged.
        // Some Segments can appear in more than one of the inEdges, outEdges,
        // inEdgesRc, outEdgesRc vectors. To avoid working with invalidated
        // edge_descriptors, we work with Segment ids instead.
        std::map<uint64_t, Segment> segmentMap;
        vector<uint64_t> inEdgeIds;
        for(const Segment e: inEdges) {
            inEdgeIds.push_back(id(e));
            segmentMap.insert(make_pair(id(e), e));
        }
        vector<uint64_t> outEdgeIds;
        for(const Segment e: outEdges) {
            outEdgeIds.push_back(id(e));
            segmentMap.insert(make_pair(id(e), e));
        }
        for(const Segment e: inEdgesRc) {
            segmentMap.insert(make_pair(id(e), e));
        }
        for(const Segment e: outEdgesRc) {
            segmentMap.insert(make_pair(id(e), e));
        }
        for(const uint64_t inEdgeId: inEdgeIds) {
            const Segment eNew = disconnectAtEnd(segmentMap.at(inEdgeId));
            segmentMap[id(eNew)] = eNew;
            const Segment eNewRc = assemblyGraph[eNew].eRc;
            segmentMap[id(eNewRc)] = eNewRc;
        }
        for(const uint64_t outEdgeId: outEdgeIds) {
            const Segment eNew = disconnectAtBeginning(segmentMap.at(outEdgeId));
            segmentMap[id(eNew)] = eNew;
            const Segment eNewRc = assemblyGraph[eNew].eRc;
            segmentMap[id(eNewRc)] = eNewRc;
        }

        // Check that all is good.
        for(const auto&[segmentId, e]: segmentMap) {
            SHASTA2_ASSERT(segmentId == id(e));
        }



        // Now we can just connect in-phase or out-of-phase inEdge/outEdge pairs.
        // Loop over the entrance/exit pairs to be connected.
        for(uint64_t i=0; i<2; i++) {
            const uint64_t entranceId = inEdgeIds[i];
            const uint64_t exitId = (isInPhase ? outEdgeIds[i] : outEdgeIds[1-i]);
            const Segment entrance = segmentMap.at(entranceId);
            const Segment exit = segmentMap.at(exitId);
            const vertex_descriptor vEntrance = target(entrance, assemblyGraph);
            const vertex_descriptor vExit = source(exit, assemblyGraph);

            // Connect them with a null edge.This will later be removed during compression.
            const auto[e, wasAdded] = add_edge(vEntrance, vExit, AssemblyGraphEdge(nextEdgeId++), assemblyGraph);
            createReverseComplementEdge(e);

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



bool AssemblyGraph::detangleAndReadFollowingIteration(const string& debugOutputBaseName)
{
    // EXPOSE WHEN CODE STABILIZES.
    const uint64_t lengthThreshold = 100000;

    performanceLog << timestamp << "AssemblyGraph::detangleAndReadFollowingIteration begins: " <<
        debugOutputBaseName << endl;
    const bool debug = false;
    if(debug) {
        cout << timestamp << "AssemblyGraph::detangleIteration begins." << endl;
    }

    AssemblyGraph& assemblyGraph = *this;

    // Map the vertices to integers.
    // This is needed below to compute connected components.
    std::map<vertex_descriptor, uint64_t> vertexIndexMap;
    std::vector<vertex_descriptor> vertexTable;
    BGL_FORALL_VERTICES(v, assemblyGraph, AssemblyGraph) {
        vertexIndexMap.insert(make_pair(v, vertexTable.size()));
        vertexTable.push_back(v);
    }

    // Compute connected components using only short edges.
    // Each non-trivial connected component will become a tangle.
    DisjointSets disjointSets(vertexTable.size());
    BGL_FORALL_EDGES(segment, assemblyGraph, AssemblyGraph) {
        if(assemblyGraph[segment].length() < lengthThreshold) {
            const vertex_descriptor v0 = source(segment, assemblyGraph);
            const vertex_descriptor v1 = target(segment, assemblyGraph);
            disjointSets.unionSet(vertexIndexMap.at(v0), vertexIndexMap.at(v1));
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
    if(debug) {
        cout << "Found " << tangles.size() << " tangles." << endl;
    }

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
    if(true) {
        ofstream csv(debugOutputBaseName + "-Tangles-Bandage.csv");
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



    // Detangle or read following on each pair of reverse complemented tangles.
    bool somethingWasDone = false;
    for(uint64_t tangleId=0; tangleId<tangles.size(); tangleId++) {
        if(tangleId <= tangleRc[tangleId]) {
            const Tangle tangle(assemblyGraph, tangles[tangleId]);

            ofstream html;
            if(debug) {

                cout << "Working on tangle " << tangleId <<
                    " and its reverse complement tangle " << tangleRc[tangleId] << endl;
                cout << "This tangle has " <<
                    tangle.entrances.size() << " entrances, " <<
                    tangle.exits.size() << " exits, " <<
                    tangle.tangleEdges.size() << " segments, " <<
                    tangle.entrances.size() << " vertices." << endl;

                // Tangle details go to the html file.
                html.open(debugOutputBaseName + "-Tangle-" + to_string(tangleId) + ".html");
                writeHtmlBegin(html, "Tangle " + to_string(tangleId));
                html << "<h1>Tangle " << tangleId << "</h1>";
                tangle.writeHtml(html);
            }

            if(tangle.entrances.empty() or tangle.exits.empty()) {
                if(debug) {
                    cout << "Skipping tangle " << tangleId << " with no entrances and/or no exits." << endl;
                }
                continue;
            }

            // Try detangling.
            bool success = detangleStrandSymmetric(tangle, html);

            // If detangling did not work, try read following.
            if(success) {
                if(debug) {
                    cout << "Detangling was successful on this tangle." << endl;
                }
            } else {
                if(debug) {
                    cout << "Detangling was not successful on this tangle. Trying read following." << endl;
                }
                success = readFollowingStrandSymmetric(tangleId, tangle, html);
                if(debug) {
                    if(success) {
                        cout << "Read following was successful on this tangle." << endl;
                    } else {
                        cout << "Read following was not successful on this tangle." << endl;
                    }
                }
            }

            if(success) {
                somethingWasDone = true;
            }
        }
    }

    performanceLog << timestamp << "AssemblyGraph::detangleAndReadFollowingIteration ends: " <<
        debugOutputBaseName << endl;

    return somethingWasDone;
}



bool AssemblyGraph::detangleStrandSymmetric(const Tangle& tangle, ostream& html)
{
    if(tangle.isSelfComplementary()) {
        if((tangle.entrances.size() == 2) and (tangle.exits.size() == 2)) {
            return detangleSelfComplementaryTangle2By2(tangle, html);
        } else {
            return detangleSelfComplementaryTangle(tangle, html);
        }
    } else {
        return detangleTanglePair(tangle, html);
    }
}



// This detangles the tangle passed in as an argument and its reverse complement.
// The tangle must no be not self-complementary,
bool AssemblyGraph::detangleTanglePair(
    const Tangle& tangle, ostream& html)
{
    AssemblyGraph& assemblyGraph = *this;
    SHASTA2_ASSERT(not tangle.isSelfComplementary());
    SHASTA2_ASSERT(not tangle.entrances.empty());
    SHASTA2_ASSERT(not tangle.exits.empty());

    if(html) {
        html << "<h2>Detangling</h2>";
    }

    // Compute the block structure of the TangleMatrix.
    vector<TangleMatrix> blocks;
    tangle.tangleMatrix().findBlockStructure(blocks);



    // Run a G-test for each block.
    // Gather entrance/exit pairs to be connected.
    vector< pair<uint64_t, uint64_t> > connectPairs;
    vector<GTest> gTests;
    uint64_t failedGTestCount = 0;
    uint64_t reliableGTestCount = 0;
    bool connectionFailure = false;
    for(const TangleMatrix& block: blocks) {
        const bool onlyConsiderPermutation = (block.entrances.size() == block.exits.size());
        const bool onlyConsiderInjective = onlyConsiderPermutation;
        const GTest& gTest = gTests.emplace_back(block.tangleMatrix, assemblyGraph.options.detangleEpsilon,
        onlyConsiderInjective, onlyConsiderPermutation);
        if(not gTest.success) {
            ++failedGTestCount;
        } else {

            // This G-test was successful.
            // Find out if it the best hypothesis is sufficiently reliable
            // accordign to our current options.
            bool isReliable = true;
            const auto& bestHypothesis = gTest.hypotheses[0];
            const double bestG = bestHypothesis.G;
            if(bestG > assemblyGraph.options.detangleMaxLogP) {
                isReliable = false;
            } else {
                if(gTest.hypotheses.size() > 1) {
                    if(gTest.hypotheses[1].G - bestG < assemblyGraph.options.detangleMinLogPDelta) {
                        // There is a second hypothesis, and it is not sufficiently well separated
                        // from the best hypothesis.
                        isReliable = false;
                    }
                }
            }
            if(isReliable) {
                ++reliableGTestCount;

                // This test is reliable. Find the pairs to be connected.
                for(uint64_t i=0; i<block.entrances.size(); i++) {
                    const Segment entrance = block.entrances[i];
                    for(uint64_t j=0; j<block.exits.size(); j++) {
                        const Segment exit = block.exits[j];
                        if(bestHypothesis.connectivityMatrix[i][j]) {
                            if(not canConnect(entrance, exit)) {
                                connectionFailure = true;
                                if(html) {
                                    html << "<br>Cannot connect " << id(entrance) <<
                                        " " << id(exit);
                                }
                            }
                            connectPairs.push_back(make_pair(id(entrance), id(exit)));
                        }
                    }
                }
            }
        }
    }



    if(html) {
        for(uint64_t iBlock=0; iBlock<blocks.size(); iBlock++) {
            const TangleMatrix& block = blocks[iBlock];
            const GTest gTest = gTests[iBlock];

            html << "<h3>Tangle matrix block " << iBlock << "</h3>";
            block.writeTotalTangleMatrix(html);

            if(gTest.success) {
                const auto& bestHypothesis = gTest.hypotheses.front();
                const double bestG = bestHypothesis.G;
                html <<
                    "<br>G-test best hypothesis:"
                    "<table><tr><th>";
                for(const AssemblyGraph::edge_descriptor exit: block.exits) {
                    html << "<th>" << id(exit);
                }
                for(uint64_t i=0; i<block.entrances.size(); i++) {
                    html << "<tr><th>" << id(block.entrances[i]);
                    for(uint64_t j=0; j<block.exits.size(); j++) {
                        html << "<td class=centered>" << int(bestHypothesis.connectivityMatrix[i][j]);
                    }
                }
                html << "</table><br>G = " << bestG;
                if(gTest.hypotheses.size() > 1) {
                    const double secondBestG = gTest.hypotheses[1].G;
                    html << "<br>Second best G = " << secondBestG;
                    html << "<br>&Delta;G = " << secondBestG - bestG;
                }

                bool isReliable = true;
                if(bestG > assemblyGraph.options.detangleMaxLogP) {
                    isReliable = false;
                } else {
                    if(gTest.hypotheses.size() > 1) {
                        if(gTest.hypotheses[1].G - bestG < assemblyGraph.options.detangleMinLogPDelta) {
                            // There is a second hypothesis, and it is not sufficiently well separated
                            // from the best hypothesis.
                            isReliable = false;
                        }
                    }
                }

                if(isReliable) {
                    html << "<br>The G-test indicates that the best hypothesis for this block is reliable.";
                } else {
                    html << "<br>The G-test indicates that the best hypothesis for this block is not reliable.";
               }

            } else {
                html << "<br>The G-test for this block failed.";
            }
        }
    }


    if(failedGTestCount > 0) {
        if(html) {
            html << "<br>Not detangling because the G-test for " << failedGTestCount <<
                " blocks failed.";
        }
        return false;
    }
    if(reliableGTestCount < blocks.size()) {
        if(html) {
            html << "<br>Not detangling because only " << reliableGTestCount <<
                " blocks out of " << blocks.size() <<
                " gave a sufficiently reliable result.";
        }
        return false;
    }
    if(connectionFailure) {
        if(html) {
            html << "<br>Not detangling because of connection failures.";
        }
        return false;
    }

    if(html) {
        html << "<br>This tangle and its reverse complement will be detangled." << endl;
    }

    // Remove all the Segments internal to this tangle
    // and their reverse complements.
    for(const Segment e: tangle.tangleEdges) {
        boost::remove_edge(assemblyGraph[e].eRc, assemblyGraph);
        boost::remove_edge(e, assemblyGraph);
    }



    // If getting here, this tangle and its reverse complement will be detangled.
    // To do this, we need to disconnect (at the beginning or end) the Segments involved.
    // To do this we use disconnectAtBeginning and disconnectAtEnd, which change
    // edge_descriptors but leave the ids unchanged.
    // Some Segments can at the same time be an entrance and/or exit
    // of this tangle or its reverse complement. To avoid working with invalidated
    // edge_descriptors, we work with Segment ids instead.
    std::map<uint64_t, Segment> segmentMap;
    vector<uint64_t> entranceIds;
    for(const Segment e: tangle.entrances) {
        entranceIds.push_back(id(e));
        segmentMap.insert(make_pair(id(e), e));
        const Segment eRc = assemblyGraph[e].eRc;
        segmentMap.insert(make_pair(id(eRc), eRc));
    }
    vector<uint64_t> exitIds;
    for(const Segment e: tangle.exits) {
        exitIds.push_back(id(e));
        segmentMap.insert(make_pair(id(e), e));
        const Segment eRc = assemblyGraph[e].eRc;
        segmentMap.insert(make_pair(id(eRc), eRc));
    }
    for(const uint64_t entranceId: entranceIds) {
        const Segment eNew = disconnectAtEnd(segmentMap.at(entranceId));
        segmentMap[id(eNew)] = eNew;
        const Segment eNewRc = assemblyGraph[eNew].eRc;
        segmentMap[id(eNewRc)] = eNewRc;
    }
    for(const uint64_t exitId: exitIds) {
        const Segment eNew = disconnectAtBeginning(segmentMap.at(exitId));
        segmentMap[id(eNew)] = eNew;
        const Segment eNewRc = assemblyGraph[eNew].eRc;
        segmentMap[id(eNewRc)] = eNewRc;
    }

    // Check that all is good.
    for(const auto&[segmentId, e]: segmentMap) {
        SHASTA2_ASSERT(segmentId == id(e));
    }


    // Make the connections.
    // This also does the same for the reverse complement of this Tangle.
    for(const auto&[entranceId, exitId]: connectPairs) {
        const AssemblyGraph::edge_descriptor entrance = segmentMap.at(entranceId);
        const AssemblyGraph::edge_descriptor exit = segmentMap.at(exitId);
        const edge_descriptor eNew = connect(entrance, exit);
        createReverseComplementEdge(eNew);
    }



    // Remove all the Tangle vertices and their reverse complements
    // that are now isolated. These are the ones that were
    // not connected to any entrance or exit.
    for(const vertex_descriptor v: tangle.tangleVertices) {
        const vertex_descriptor vRc = assemblyGraph[v].vRc;

        if((in_degree(v, assemblyGraph) == 0) and (out_degree(v, assemblyGraph) == 0)) {
            boost::remove_vertex(v, assemblyGraph);
        }

        if((in_degree(vRc, assemblyGraph) == 0) and (out_degree(vRc, assemblyGraph) == 0)) {
            boost::remove_vertex(vRc, assemblyGraph);
        }
    }

    return true;
}



bool AssemblyGraph::detangleSelfComplementaryTangle(
    const Tangle& tangle,
    ostream& html)
{
    SHASTA2_ASSERT(tangle.isSelfComplementary());

    if(html) {
        html << "<br>This tangle will not be detangled because "
            "detangleSelfComplementaryTangle is not implemented." << endl;
    }

    return false;
}



bool AssemblyGraph::detangleSelfComplementaryTangle2By2(
    const Tangle& tangle,
    ostream& html)
{
    AssemblyGraph& assemblyGraph = *this;
    if(html) {
        html << "<br>AssemblyGraph::detangleSelfComplementaryTangle2By2 begins" << endl;
    }

    const vector<Segment>& entrances = tangle.entrances;
    const vector<Segment>& exits = tangle.exits;

    SHASTA2_ASSERT(tangle.isSelfComplementary());
    SHASTA2_ASSERT(entrances.size() == 2);
    SHASTA2_ASSERT(exits.size() == 2);

    // In a self-complementary tangle, each exit is the reverse complement
    // of an entrance. We can't connect an entrance with its own reverse
    // complement, so each entrance can only possibly be connected to
    // the one exit with is not its reverse complement.
    // Therefore we don't need to compute a TangleMatrix and run a G-test.
    // We just need to check that that connection is possible.

    // Store the reverse complement of each entrance.
    vector<Segment> entrancesRc(2);
    for(uint64_t i=0; i<2; i++) {
        entrancesRc[i] = assemblyGraph[entrances[i]].eRc;
        SHASTA2_ASSERT(std::ranges::contains(exits, entrancesRc[i]));
    }

    // Define the entrance/exit pair that we will explictly connect.
    // The remaining entrance/exit pair will be connected
    // automatically.
    const Segment entrance = entrances[0];
    const Segment exit = entrancesRc[1];
    const uint64_t entranceId = id(entrance);
    const uint64_t exitId = id(exit);
    if(html) {
        html << "<br>Will attempt to connect " << entranceId << " with " << exitId << endl;
    }

    // Check if  we can connect this entrance/exit pair.
    if(not assemblyGraph.canConnect(entrance, exit)) {
        if(html) {
            html << "<br>Not detangling because can't connect " << entranceId <<
                " with " << exitId << endl;
        }
        return false;
    }

    if(html) {
        html << "<br>This tangle will be detangled." << endl;
    }

    // Remove all the Segments internal to this tangle.
    for(const Segment e: tangle.tangleEdges) {
        boost::remove_edge(e, assemblyGraph);
    }



    // If getting here, this tangle will be detangled.
    // To do this, we need to disconnect (at the beginning or end) the Segments involved.
    // To do this we use disconnectAtBeginning and disconnectAtEnd, which change
    // edge_descriptors but leave the ids unchanged.
    // Some Segments can at the same time be an entrance and/or exit
    // of this tangle or its reverse complement. To avoid working with invalidated
    // edge_descriptors, we work with Segment ids instead.
    std::map<uint64_t, Segment> segmentMap;
    vector<uint64_t> entranceIds;
    for(const Segment e: tangle.entrances) {
        entranceIds.push_back(id(e));
        segmentMap.insert(make_pair(id(e), e));
        const Segment eRc = assemblyGraph[e].eRc;
        segmentMap.insert(make_pair(id(eRc), eRc));
    }
    vector<uint64_t> exitIds;
    for(const Segment e: tangle.exits) {
        exitIds.push_back(id(e));
        segmentMap.insert(make_pair(id(e), e));
        const Segment eRc = assemblyGraph[e].eRc;
        segmentMap.insert(make_pair(id(eRc), eRc));
    }
    for(const uint64_t entranceId: entranceIds) {
        const Segment eNew = disconnectAtEnd(segmentMap.at(entranceId));
        segmentMap[id(eNew)] = eNew;
        const Segment eNewRc = assemblyGraph[eNew].eRc;
        segmentMap[id(eNewRc)] = eNewRc;
    }
    for(const uint64_t exitId: exitIds) {
        const Segment eNew = disconnectAtBeginning(segmentMap.at(exitId));
        segmentMap[id(eNew)] = eNew;
        const Segment eNewRc = assemblyGraph[eNew].eRc;
        segmentMap[id(eNewRc)] = eNewRc;
    }

    // Check that all is good.
    for(const auto&[segmentId, e]: segmentMap) {
        SHASTA2_ASSERT(segmentId == id(e));
    }

    // Make the connection and its reverse complement..
    const edge_descriptor eNew = connect(segmentMap.at(entranceId), segmentMap.at(exitId));
    createReverseComplementEdge(eNew);

    // Remove all the Tangle vertices that are now isolated.
    // These are the ones that were not connected to any entrance or exit.
    for(const vertex_descriptor v: tangle.tangleVertices) {
        if((in_degree(v, assemblyGraph) == 0) and (out_degree(v, assemblyGraph) == 0)) {
            boost::remove_vertex(v, assemblyGraph);
        }
    }

    return true;
}


