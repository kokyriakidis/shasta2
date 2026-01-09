
uint64_t AssemblyGraph::prune()
{
    AssemblyGraph& assemblyGraph = *this;
    const bool debug = false;

    uint64_t totalPruneCount = 0;

    vector<vertex_descriptor> verticesToBeRemoved;
    vector<edge_descriptor> edgesToBeRemoved;
    for(uint64_t iteration=0; iteration<options.pruneIterationCount; iteration++) {

        // Edge pruning for this iteration.
        edgesToBeRemoved.clear();
        uint64_t pruneCount = 0;
        uint64_t prunedLength = 0;
        BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
            const uint64_t offset = assemblyGraph[e].offset();

            // If long enough, don't prune it.
            if(offset > options.pruneLength) {
                continue;
            }

            // See if it can be pruned.
            const vertex_descriptor v0 = source(e, assemblyGraph);
            const vertex_descriptor v1 = target(e, assemblyGraph);
            if((in_degree(v0, assemblyGraph) == 0) or (out_degree(v1, assemblyGraph) == 0)) {
                edgesToBeRemoved.push_back(e);
                ++pruneCount;
                prunedLength += offset;
            }
        }

        // Remove the edges.
        for(const edge_descriptor e: edgesToBeRemoved) {
            if(debug) {
                cout << "Pruned " << assemblyGraph[e].id << ", offset " << assemblyGraph[e].offset() << endl;
            }
            boost::remove_edge(e, assemblyGraph);
        }

        // Now remove any vertices that are left isolated.
        verticesToBeRemoved.clear();
        BGL_FORALL_VERTICES(v, assemblyGraph, AssemblyGraph) {
            if((in_degree(v, assemblyGraph) == 0) and (out_degree(v, assemblyGraph) == 0)) {
                verticesToBeRemoved.push_back(v);
            }
        }
        for(const vertex_descriptor v: verticesToBeRemoved) {
            boost::remove_vertex(v, assemblyGraph);
        }

        cout << "Prune iteration removed " << pruneCount <<
            " edges with total estimated length " << prunedLength << endl;

        if(pruneCount == 0) {
            break;
        }

        totalPruneCount += pruneCount;
    }

    return totalPruneCount;
}




uint64_t AssemblyGraph::detangleVertices(Detangler& detangler)
{
    uint64_t changeCount = 0;

    for(uint64_t iteration=0; iteration<options.detangleMaxIterationCount; iteration++) {
        const uint64_t iterationChangeCount = detangleVerticesIteration(detangler);
        if(iterationChangeCount > 0) {
            changeCount += iterationChangeCount;
            changeCount += compress();
            cout << "Detangle vertices iteration " << iteration << ": " << changeCount <<
                " successful detangling operations." << endl;
        } else {
            break;
        }
    }

    return changeCount;
}



uint64_t AssemblyGraph::detangleVerticesIteration(Detangler& detangler)
{
    AssemblyGraph& assemblyGraph = *this;

    // Gather vertices on which we will attempt detangling.
    // Each generates a tangle with just one vertex.
    vector< vector<vertex_descriptor> > detanglingCandidates;
    BGL_FORALL_VERTICES(v, assemblyGraph, AssemblyGraph) {

        if(
            (in_degree(v, assemblyGraph) > 1) and
            (out_degree(v, assemblyGraph) > 1)
             ) {
            detanglingCandidates.emplace_back(vector<vertex_descriptor>({v}));
        }

    }

    // Do the detangling.
    return detangleLowLevel(detanglingCandidates, detangler);
}



uint64_t AssemblyGraph::detangleEdges(Detangler& detangler)
{
    uint64_t changeCount = 0;

    for(uint64_t iteration=0; iteration<options.detangleMaxIterationCount; iteration++) {
        const uint64_t iterationChangeCount = detangleEdgesIteration(detangler);
        if(iterationChangeCount > 0) {
            changeCount += iterationChangeCount;
            changeCount += compress();
            cout << "Detangle edges iteration " << iteration << ": " << iterationChangeCount <<
                " successful detangling operations." << endl;
        } else {
            break;
        }
    }

    return changeCount;
}



uint64_t AssemblyGraph::detangleEdgesIteration(Detangler& detangler)
{
    AssemblyGraph& assemblyGraph = *this;
    // cout << "Edge detangling begins." << endl;

    // Gather edges on which we will attempt detangling.
    // Each generates a tangle with just two vertices.
    vector< vector<vertex_descriptor> > detanglingCandidates;
    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        const vertex_descriptor v0 = source(e, assemblyGraph);
        const vertex_descriptor v1 = target(e, assemblyGraph);

        const uint64_t in0 = in_degree(v0, assemblyGraph);
        const uint64_t out0 = out_degree(v0, assemblyGraph);
        const uint64_t in1 = in_degree(v1, assemblyGraph);
        const uint64_t out1 = out_degree(v1, assemblyGraph);

        const bool isTangleEdge = (in0 > 1) and (out0 == 1) and (in1 == 1) and (out1 > 1);

        // For cross edges, don't use in-degree and out-degree.
        // This is needed to make sure that bubble edges are not classified as cross-edges.
        const bool isCrossEdge =
            (countDistinctTargetVertices(v0) > 1) and
            (countDistinctSourceVertices(v1) > 1);

        const bool isShort = assemblyGraph[e].offset() <= options.detangleMaxCrossEdgeLength;

        // Tangle edges without length limitations.
        // Cross edges limited by length.
        if(isTangleEdge or(isShort and isCrossEdge)) {
            detanglingCandidates.emplace_back(vector<vertex_descriptor>({v0, v1}));
        }
    }

    // Do the detangling.
    return detangleLowLevel(detanglingCandidates, detangler);
}



uint64_t AssemblyGraph::detangleLowLevel(
    const vector< vector<vertex_descriptor> >& detanglingCandidates,
    Detangler& detangler)
{
    const bool debug = false;

    AssemblyGraph& assemblyGraph = *this;

    std::set<vertex_descriptor> removedVertices;
    uint64_t attemptCount = 0;
    uint64_t successCount = 0;
    for(const vector<vertex_descriptor>& tangleVertices: detanglingCandidates) {

        // If any of the vertices in this tangle have been removed, by previous
        // detangling operations, skip it.
        bool skip = false;
        for(const vertex_descriptor v: tangleVertices) {
            if(removedVertices.contains(v)) {
                skip = true;
                break;
            }
        }
        if(skip) {
            continue;
        }



        // Attempt detangling for the tangle defined by these vertices.
        ++attemptCount;
        Tangle1 tangle(assemblyGraph, tangleVertices);
        if(debug) {
            const TangleMatrix1& tangleMatrix = tangle.tangleMatrix();
            cout << "Tangle with " << tangleMatrix.entrances.size() << " entrances and " <<
                tangleMatrix.exits.size() << " exits." << endl;

            cout << "Entrances:";
            for(const edge_descriptor e: tangleMatrix.entrances) {
                cout << " " << assemblyGraph[e].id;
            }
            cout << endl;

            cout << "Exits:";
            for(const edge_descriptor e: tangleMatrix.exits) {
                cout << " " << assemblyGraph[e].id;
            }
            cout << endl;

            for(uint64_t iEntrance=0; iEntrance<tangleMatrix.entrances.size(); iEntrance++) {
                for(uint64_t iExit=0; iExit<tangleMatrix.exits.size(); iExit++) {
                    cout << tangleMatrix.tangleMatrix[iEntrance][iExit] << " ";
                }
            }
            cout << endl;
        }

        const bool success = detangler(tangle);
        if(success) {
            if(debug) {
                cout << "Detangle was successful." << endl;
            }
            for(vertex_descriptor v: tangle.removedVertices) {
                removedVertices.insert(v);
            }
            ++successCount;
        } else {
            if(debug) {
                cout << "Detangle failed." << endl;
            }
        }


    }

    // The detangling process can generate empty edges. Remove them.
    removeEmptyEdges();

    // cout << "Attempted detangling for " << attemptCount << " tangles." << endl;
    // cout << "Detangling was successful for " << successCount << " tangles." << endl;
    return successCount;
}



uint64_t AssemblyGraph::detangleHighLevel(Detangler& detangler)
{
    performanceLog << timestamp << "AssemblyGraph::detangle begins." << endl;

    const uint64_t verticesChangeCount = detangleVertices(detangler);
    const uint64_t edgesChangeCount = detangleEdges(detangler);
    const uint64_t shortTanglesChangeCount = detangleShortTangles(detangler);

    const uint64_t changeCount = verticesChangeCount + edgesChangeCount + shortTanglesChangeCount;

    performanceLog << timestamp << "AssemblyGraph::detangle ends." << endl;
    return changeCount;
}




// Detangle short tangles.
uint64_t AssemblyGraph::detangleShortTangles(Detangler& detangler)
{
    const AssemblyGraph& assemblyGraph = *this;

    // EXPOSE WHEN CODE STABILIZES.
    const uint64_t shortTangleMaxStepCount = 10;

    // Map vertices to integers.
    std::map<vertex_descriptor, uint64_t> vertexIndexMap;
    vector<vertex_descriptor> vertexTable;
    uint64_t vertexIndex = 0;
    BGL_FORALL_VERTICES(v, assemblyGraph, AssemblyGraph) {
        vertexIndexMap.insert(make_pair(v, vertexIndex++));
        vertexTable.push_back(v);
    }
    const uint64_t n = vertexIndexMap.size();

    // Initialize the disjoint set data structures.
    vector<uint64_t> rank(n);
    vector<uint64_t> parent(n);
    boost::disjoint_sets<uint64_t*, uint64_t*> disjointSets(&rank[0], &parent[0]);
    for(uint64_t i=0; i<n; i++) {
        disjointSets.make_set(i);
    }



    // Compute connected components of the AssemblyGraph using only short edges.
    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {

        // If not a short edge, skip it.
        if(assemblyGraph[e].size() > shortTangleMaxStepCount) {
            continue;
        }

        // Find the index corresponding to the source vertex.
        const vertex_descriptor v0 = source(e, assemblyGraph);
        const auto it0 = vertexIndexMap.find(v0);
        SHASTA2_ASSERT(it0 != vertexIndexMap.end());
        const uint64_t i0 = it0->second;

        // Find the index corresponding to the target vertex.
        const vertex_descriptor v1 = target(e, assemblyGraph);
        const auto it1 = vertexIndexMap.find(v1);
        SHASTA2_ASSERT(it1 != vertexIndexMap.end());
        const uint64_t i1 = it1->second;

        // Update the disjoint sets for this edge.
        disjointSets.union_set(i0, i1);
    }



    // Gather the vertices in each connected component.
    vector< vector<uint64_t> > components(n);
    for(uint64_t i=0; i<n; i++) {
        components[disjointSets.find_set(i)].push_back(i);
    }



    // Each connected components with size greater than 1 generates a short tangle.
    vector< vector<vertex_descriptor> > detanglingCandidates;
    for(const vector<uint64_t>& component: components) {
        if(component.size() <= 1) {
            continue;
        }

        // Gather the tangle vertices.
        vector<vertex_descriptor> tangleVertices;
        for(const uint64_t i: component) {
            tangleVertices.push_back(vertexTable[i]);
        }

        // Add it to our detangling candidates.
        detanglingCandidates.push_back(tangleVertices);
    }

    // Do the detangling.
    const uint64_t successCount = detangleLowLevel(detanglingCandidates, detangler);

    cout << "Of " << detanglingCandidates.size() <<
        " short tangles, " << successCount << " were detangled successfully." << endl;
    compress();

    return successCount;
}






// Simplify Superbubbles by turning them into bubbles via clustering
// of oriented read journeys.
uint64_t AssemblyGraph::simplifySuperbubbles(Detangler& detangler)
{

    // Find the superbubbles, then remove superbubbles that are entirely
    // contained in another superbubble.
    vector<Superbubble> superbubbles;
    findSuperbubbles(superbubbles);
    removeContainedSuperbubbles(superbubbles);

    // Count the number of true Superbubbles, excluding bubbles.
    uint64_t count = 0;
    for(const Superbubble& superbubble: superbubbles) {
        if(not superbubble.isBubble()) {
            ++count;
        }
    }
    cout << "Found " << superbubbles.size() << " superbubbles of which " <<
        count << " are not simple bubbles." << endl;

    uint64_t simplifiedCount = 0;
    for(const Superbubble& superbubble: superbubbles) {
        if(not superbubble.isBubble()) {

            // First try simplify by detangling.
            bool success = simplifySuperbubbleByDetangling(superbubble, detangler);

            // If that did not work, try simplify by clustering.
            if(not success) {
                simplifySuperbubbleByClustering(superbubble,
                    options.simplifySuperbubbleMinCoverage,
                    options.simplifySuperbubbleMaxOffset);
            }

            // If one of the above was successful, increment our counter.
            if(success) {
                ++simplifiedCount;
            }
        }
    }

    // The detangling process can generate empty edges. Remove them.
    removeEmptyEdges();

    return simplifiedCount;
}





bool AssemblyGraph::simplifySuperbubbleByDetangling(
    const Superbubble& superbubble,
    Detangler& detangler)
{
    AssemblyGraph& assemblyGraph = *this;

    const bool debug = false;

    Tangle1 tangle(assemblyGraph, superbubble.internalVertices);

    if(debug) {
        cout << "simplifySuperbubbleByDetangling begins." << endl;

        const TangleMatrix1& tangleMatrix = tangle.tangleMatrix();
        cout << "Tangle with " << tangleMatrix.entrances.size() << " entrances and " <<
            tangleMatrix.exits.size() << " exits." << endl;

        cout << "Entrances:";
        for(const edge_descriptor e: tangleMatrix.entrances) {
            cout << " " << assemblyGraph[e].id;
        }
        cout << endl;

        cout << "Exits:";
        for(const edge_descriptor e: tangleMatrix.exits) {
            cout << " " << assemblyGraph[e].id;
        }
        cout << endl;

        for(uint64_t iEntrance=0; iEntrance<tangleMatrix.entrances.size(); iEntrance++) {
            for(uint64_t iExit=0; iExit<tangleMatrix.exits.size(); iExit++) {
                cout << tangleMatrix.tangleMatrix[iEntrance][iExit] << " ";
            }
        }
        cout << endl;
    }

    const bool success = detangler(tangle);
    if(success) {
        if(debug) {
            cout << "Detangle was successful." << endl;
        }
    } else {
        if(debug) {
            cout << "Detangle failed." << endl;
        }
    }

    return success;

}



bool AssemblyGraph::simplifySuperbubbleByClustering(
    const Superbubble& superbubble,
    uint64_t minCoverage,
    uint64_t maxOffset)
{
    AssemblyGraph& assemblyGraph = *this;

    const AnchorId anchorIdA = assemblyGraph[superbubble.sourceVertex].anchorId;
    const AnchorId anchorIdB = assemblyGraph[superbubble.targetVertex].anchorId;

    const bool debug = false;



    if(debug) {
        cout << "Working on a superbubble with source vertex " <<
            assemblyGraph[superbubble.sourceVertex].id <<
            " and target vertex " << assemblyGraph[superbubble.targetVertex].id << endl;

        cout << "Source edges:";
        for(const edge_descriptor e: superbubble.sourceEdges) {
            cout << " " << assemblyGraph[e].id << " ";
        }
        cout << endl;

        cout << "Target edges:";
        for(const edge_descriptor e: superbubble.targetEdges) {
            cout << " " << assemblyGraph[e].id;
        }
        cout << endl;

        cout << "This superbubble consists of the following " <<
            superbubble.internalEdges.size() << " edges:" << endl;
        for(const edge_descriptor e: superbubble.internalEdges) {
            cout << assemblyGraph[e].id << " ";
        }
        cout << endl;
    }



    // Create an AnchorPair between anchorIdA and anchorIdB
    // using all oriented reads in common between anchorIdA and anchorIdB.
    AnchorPair anchorPair(anchors, anchorIdA, anchorIdB, false);

    if(debug) {
        cout << "The initial anchor pair " <<
            anchorIdToString(anchorIdA) << " " <<
            anchorIdToString(anchorIdB) <<
            " has " << anchorPair.orientedReadIds.size() <<
            " oriented reads." << endl;
    }

    // We only want to use OrientedReadIds that appear at least once in the internal
    // edges of the superbubble.
    vector<OrientedReadId> allowedOrientedReadIds;
    for(const edge_descriptor e: superbubble.internalEdges) {
        const AssemblyGraphEdge& edge = assemblyGraph[e];
        for(const AssemblyGraphEdgeStep& step: edge) {
            const AnchorPair& stepAnchorPair = step.anchorPair;
            copy(stepAnchorPair.orientedReadIds.begin(), stepAnchorPair.orientedReadIds.end(),
                back_inserter(allowedOrientedReadIds));
        }
    }
    deduplicate(allowedOrientedReadIds);
    if(debug) {
        cout << "The internal edges of the superbubble contain " <<
            allowedOrientedReadIds.size() << " oriented reads." << endl;
    }

    // Only keep OrientedReadIds in the allowed set.
    vector<OrientedReadId> newOrientedReadIds;
    std::set_intersection(
        anchorPair.orientedReadIds.begin(), anchorPair.orientedReadIds.end(),
        allowedOrientedReadIds.begin(), allowedOrientedReadIds.end(),
        back_inserter(newOrientedReadIds));
    anchorPair.orientedReadIds.swap(newOrientedReadIds);

    if(debug) {
        cout << "The final anchor pair has " << anchorPair.orientedReadIds.size() <<
            " oriented reads." << endl;
    }

    // Cluster the oriented reads in the AnchorPair.
    vector<AnchorPair> newAnchorPairs;
    anchorPair.splitByClustering(anchors, journeys, options.clusteringMinJaccard, newAnchorPairs);


    if(debug) {
        cout << "Found " << newAnchorPairs.size() << " split AnchorPairs:" << endl;
        for(const AnchorPair& newAnchorPair: newAnchorPairs) {
            cout << "AnchorPair with coverage " << newAnchorPair.orientedReadIds.size() <<
                ", offset " << newAnchorPair.getAverageOffset(anchors) << endl;
        }

#if 0
        // Also write out assembled sequences.
        for(uint64_t i=0; i<newAnchorPairs.size(); i++) {
            const AnchorPair& newAnchorPair = newAnchorPairs[i];
            ostream html(0);
            LocalAssembly2 localAssembly(
                anchors, html, false,
                options.aDrift,
                options.bDrift,
                newAnchorPair);
            localAssembly.run(false, options.maxAbpoaLength);
            vector<shasta2::Base> sequence;
            localAssembly.getSequence(sequence);
            cout << ">" << i << endl;
            copy(sequence.begin(), sequence.end(), ostream_iterator<shasta2::Base>(cout));
            cout << endl;
        }
#endif
    }

    // Only keep the ones with coverage at least minCoverage.
    for(uint64_t i=0; i<newAnchorPairs.size(); i++) {
        if(newAnchorPairs[i].orientedReadIds.size() < minCoverage) {
            newAnchorPairs.resize(i);
            break;
        }
    }

    if(debug) {
        cout << "Kept " << newAnchorPairs.size() << " split AnchorPairs with sizes:";
        for(const AnchorPair& newAnchorPair: newAnchorPairs) {
            cout << " " << newAnchorPair.orientedReadIds.size();
        }
        cout << endl;
    }

    // If there are no AnchorPairs with sufficient coverage, we can't simplify this Superbubble.
    if(newAnchorPairs.empty()) {
        if(debug) {
            cout << "This superbubble cannot be simplified because there are no usable anchor pairs." << endl;
        }
        return false;
    }

    // If any of these AnchorPairs have a long offset, don't do it.
    // This only works well for small superbubbles, and for long superbubbles
    // it can destroy correct sequence.
    for(const AnchorPair& newAnchorPair: newAnchorPairs) {
        if(newAnchorPair.getAverageOffset(anchors) > maxOffset) {
            if(debug) {
                cout << "Skipping this superbubble due to large offset." << endl;
            }
            return false;
        }
    }

    // We replace the Superbubble with a bubble created using these new AnchorPairs.
    // Each of the new AnchorPairs we kept will generate a branch of the new bubble.
    for(const AnchorPair& newAnchorPair: newAnchorPairs) {
        const uint64_t offset = newAnchorPair.getAverageOffset(anchors);
        edge_descriptor e;
        bool edgeWasAdded;
        tie(e, edgeWasAdded) = add_edge(
            superbubble.sourceVertex, superbubble.targetVertex,
            AssemblyGraphEdge(nextEdgeId++), assemblyGraph);
        AssemblyGraphEdge& edge = assemblyGraph[e];
        edge.emplace_back(newAnchorPair, offset);
        if(debug) {
            cout << "Created new edge " << edge.id << endl;
        }
    }

    // Now we can remove all the internal edges of the Superbubble.
    for(const edge_descriptor e: superbubble.internalEdges) {
        if(debug) {
            cout << "Removed edge " << assemblyGraph[e].id << endl;
        }
        boost::remove_edge(e, assemblyGraph);
    }

    // Also remove the internal vertices of the Superbubble.
    for(const vertex_descriptor v: superbubble.internalVertices) {
        SHASTA2_ASSERT(in_degree(v, assemblyGraph) == 0);
        SHASTA2_ASSERT(out_degree(v, assemblyGraph) == 0);
        boost::remove_vertex(v, assemblyGraph);
    }

    return true;
}





// Local search that continues as long as we have only one way to move.
void AssemblyGraph::forwardLocalSearch(
    edge_descriptor eStart,
    uint64_t lowCoverageThreshold,
    uint64_t highCoverageThreshold,
    vector<edge_descriptor>& edges) const
{
    const bool debug = false;
    const AssemblyGraph& assemblyGraph = *this;

    edges.clear();
    edge_descriptor e0 = eStart;

    vector<edge_descriptor> entrances(1, eStart);
    vector<edge_descriptor> exits;

    // Main iteration.
    while(true) {
        const vertex_descriptor v0 = target(e0, assemblyGraph);
        exits.clear();
        BGL_FORALL_OUTEDGES(v0, e1, assemblyGraph, AssemblyGraph) {
            exits.push_back(e1);
        }

        TangleMatrix tangleMatrix(assemblyGraph, entrances, exits, 0, options.aDrift, options.bDrift);

        if(debug) {
            cout << "Starting from " << assemblyGraph[e0].id << " found:" << endl;
            for(uint64_t iExit=0; iExit<exits.size(); iExit++) {
                cout << assemblyGraph[exits[iExit]].id << " " <<
                    tangleMatrix.tangleMatrix[0][iExit].orientedReadIds.size() << endl;
            }
        }

        // Counts exits by type.
        // Significant: coverage >= highCoverageThreshold
        // Insignificant: coverage <= lowCoverageThreshold
        // Ambiguous: lowCoverageThreshold< coverage < highCoverageThreshold.
        uint64_t significantCount = 0;
        uint64_t insignificantCount = 0;
        uint64_t ambiguousCount = 0;
        for(uint64_t iExit=0; iExit<exits.size(); iExit++) {
            const uint64_t coverage = tangleMatrix.tangleMatrix[0][iExit].orientedReadIds.size();
            if(coverage <= lowCoverageThreshold) {
                ++insignificantCount;
            } else if(coverage >= highCoverageThreshold) {
                ++significantCount;
            } else {
                ++ambiguousCount;
            }
        }

        // Only keep going if we have exactly one significant exit and no ambiguous exits.
        if(significantCount != 1) {
            return;
        }
        if(ambiguousCount > 0) {
            return;
        }

        // Find the one and only significant exit.
        edge_descriptor e1;
        for(uint64_t iExit=0; iExit<exits.size(); iExit++) {
            if(tangleMatrix.tangleMatrix[0][iExit].orientedReadIds.size() >= highCoverageThreshold) {
                e1 = exits[iExit];
                break;
            }
        }

        // Check for loops.
        if(e1 == eStart) {
            break;
        }
        bool loopDetected = false;
        for(const edge_descriptor e: edges) {
            if(e == e1) {
                loopDetected = true;
                break;
            }
        }
        if(loopDetected) {
            break;
        }

        // Add it to our output edges and continue from here.
        edges.push_back(e1);
        e0 = e1;
    }
}



void AssemblyGraph::backwardLocalSearch(
    edge_descriptor eStart,
    uint64_t lowCoverageThreshold,
    uint64_t highCoverageThreshold,
    vector<edge_descriptor>& edges) const
{
    const bool debug = false;
    const AssemblyGraph& assemblyGraph = *this;

    edges.clear();
    edge_descriptor e0 = eStart;

    vector<edge_descriptor> entrances;
    vector<edge_descriptor> exits(1, eStart);

    // Main iteration.
    while(true) {
        const vertex_descriptor v0 = source(e0, assemblyGraph);
        entrances.clear();
        BGL_FORALL_INEDGES(v0, e1, assemblyGraph, AssemblyGraph) {
            entrances.push_back(e1);
        }

        TangleMatrix tangleMatrix(assemblyGraph, entrances, exits, 0, options.aDrift, options.bDrift);

        if(debug) {
            cout << "Starting from " << assemblyGraph[e0].id << " found:" << endl;
            for(uint64_t iEntrance=0; iEntrance<entrances.size(); iEntrance++) {
                cout << assemblyGraph[entrances[iEntrance]].id << " " <<
                    tangleMatrix.tangleMatrix[iEntrance][0].orientedReadIds.size() << endl;
            }
        }

        // Counts exits by type.
        // Significant: coverage >= highCoverageThreshold
        // Insignificant: coverage <= lowCoverageThreshold
        // Ambiguous: lowCoverageThreshold< coverage < highCoverageThreshold.
        uint64_t significantCount = 0;
        uint64_t insignificantCount = 0;
        uint64_t ambiguousCount = 0;
        for(uint64_t iEntrance=0; iEntrance<entrances.size(); iEntrance++) {
            const uint64_t coverage = tangleMatrix.tangleMatrix[iEntrance][0].orientedReadIds.size();
            if(coverage <= lowCoverageThreshold) {
                ++insignificantCount;
            } else if(coverage >= highCoverageThreshold) {
                ++significantCount;
            } else {
                ++ambiguousCount;
            }
        }

        // Only keep going if we have exactly one significant exit and no ambiguous exits.
        if(significantCount != 1) {
            return;
        }
        if(ambiguousCount > 0) {
            return;
        }

        // Find the one and only good entrance.
        edge_descriptor e1;
        for(uint64_t iEntrance=0; iEntrance<entrances.size(); iEntrance++) {
            if(tangleMatrix.tangleMatrix[iEntrance][0].orientedReadIds.size() >= highCoverageThreshold) {
                e1 = entrances[iEntrance];
                break;
            }
        }

        // Check for loops.
        if(e1 == eStart) {
            break;
        }
        bool loopDetected = false;
        for(const edge_descriptor e: edges) {
            if(e == e1) {
                loopDetected = true;
                break;
            }
        }
        if(loopDetected) {
            break;
        }

        // Add it to our output edges and continue from here.
        edges.push_back(e1);
        e0 = e1;
    }
}



void AssemblyGraph::testLocalSearch(
    uint64_t edgeIdStart,
    uint64_t direction,
    uint64_t lowCoverageThreshold,
    uint64_t highCoverageThreshold) const
{
    const AssemblyGraph& assemblyGraph = *this;

    edge_descriptor eStart;
    bool found = false;
    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        if(assemblyGraph[e].id == edgeIdStart) {
            eStart = e;
            found = true;
            break;
        }
    }
    if(not found) {
        cout << "Edge with id " << edgeIdStart << " does not exist." << endl;
        return;
    }

    vector<edge_descriptor> edges;
    if(direction == 0) {
        forwardLocalSearch(eStart, lowCoverageThreshold, highCoverageThreshold, edges);
    } else {
        backwardLocalSearch(eStart, lowCoverageThreshold, highCoverageThreshold, edges);
    }

    cout << "Found " << edges.size() << " edges:" << endl;
    for(const edge_descriptor e: edges) {
        cout << assemblyGraph[e].id << " ";
    }
    cout << endl;


}



void AssemblyGraph::createSearchGraph(
    uint64_t lowCoverageThreshold,
    uint64_t highCoverageThreshold)
{
    using shasta2::SearchGraph;
    AssemblyGraph& assemblyGraph = *this;

    // Create the SearchGraph.
    SearchGraph searchGraph(*this, lowCoverageThreshold, highCoverageThreshold);

    // Compute connected components.
    vector<SearchGraph> components;
    searchGraph.computeConnectedComponents(components);
    // cout << "Found " << components.size() << " non-trivial connected components of the SearchGraph." << endl;



    // Process each connected component separately.
    ofstream csv("SearchGraph-Chains.csv");
    ofstream csvBandage("SearchGraph-Bandage.csv");
    csvBandage << "Segment,Color\n";
    vector< vector<vertex_descriptor> > chains;
    for(uint64_t componentId=0; componentId<components.size(); componentId++) {
        // cout << "Working on component " << componentId << " of " << components.size() << endl;
        SearchGraph& component = components[componentId];
        transitiveReductionAny(component);
        component.removeBranches();
        // component.writeGraphviz("SearchGraph-" + to_string(componentId) + ".dot");

        findLinearVertexChains(component, chains);

        // For each chain, generate a new AssemblyGraphEdge,
        // without connecting it to the rest of the AssemblyGraph for now.
        for(const vector<vertex_descriptor>& chain: chains) {
            if(chain.size() < 2) {
                continue;
            }

            // Get some information about the first and last AssemblyGraph edge in the chain.
            const SearchGraph::vertex_descriptor sv0 = chain.front();
            const SearchGraph::vertex_descriptor sv1 = chain.back();
            const edge_descriptor e0 = component[sv0].e;
            const edge_descriptor e1 = component[sv1].e;

            const vertex_descriptor v0 = source(e0, assemblyGraph);
            const vertex_descriptor v1 = target(e1, assemblyGraph);

#if 0
            // Create new vertices for the new AssemblyGraphEdge, so it
            // will stay isolated from the rest of the AssemblyGraph for now.
            const AssemblyGraphEdge& edge0 = assemblyGraph[e0];
            const AssemblyGraphEdge& edge1 = assemblyGraph[e1];
            const AnchorId anchorId0 = edge0.firstAnchorId();
            const AnchorId anchorId1 = edge1.lastAnchorId();
            const vertex_descriptor v0 = add_vertex(AssemblyGraphVertex(anchorId0, nextVertexId++), assemblyGraph);
            const vertex_descriptor v1 = add_vertex(AssemblyGraphVertex(anchorId1, nextVertexId++), assemblyGraph);
#endif

            // Create the new AssemblyGraphEdge.
            edge_descriptor eNew;
            tie(eNew, ignore) = add_edge(v0, v1, AssemblyGraphEdge(nextEdgeId++), assemblyGraph);
            AssemblyGraphEdge& edgeNew = assemblyGraph[eNew];
            csvBandage << edgeNew.id << ",Red\n";



            // Concatenate the AssemblyGraphEdges of the chain, adding bridge steps where needed.
            for(uint64_t i=0; i<chain.size(); i++) {
                const SearchGraph::vertex_descriptor v = chain[i];
                const edge_descriptor e = component[v].e;
                const AssemblyGraphEdge& edge = assemblyGraph[e];

                // If necessary, add an AssemblyGraphEdgeStep to bridge.
                if(i != 0) {
                    const SearchGraph::vertex_descriptor vPrevious = chain[i - 1];
                    const edge_descriptor ePrevious = component[vPrevious].e;
                    const AssemblyGraphEdge& edgePrevious = assemblyGraph[ePrevious];
                    const AnchorId anchorIdPrevious = edgePrevious.lastAnchorId();
                    const AnchorId anchorIdNext = edge.firstAnchorId();
                    if(anchorIdPrevious != anchorIdNext) {
                        const AnchorPair bridgeAnchorPair = anchors.bridge(
                            edgePrevious.back().anchorPair,
                            edge.front().anchorPair,
                            options.aDrift, options.bDrift);
                        const uint64_t offset = bridgeAnchorPair.getAverageOffset(anchors);
                        edgeNew.push_back(AssemblyGraphEdgeStep(bridgeAnchorPair, offset));
                    }
                }

                // Now we can append this edge to the new edge.
                copy(edge.begin(), edge.end(), back_inserter(edgeNew));

            }

            // Now remove the edges of the chain.
            for(const SearchGraph::vertex_descriptor sv: chain) {
                const edge_descriptor e = component[sv].e;
                boost::remove_edge(e, assemblyGraph);
            }


            // Write this chain to the csv file.
            csv << edgeNew.id << ",";
            for(const vertex_descriptor v: chain) {
                const AssemblyGraph::edge_descriptor e = searchGraph[v].e;
                csv << assemblyGraph[e].id << ",";
            }
            csv << "\n";
        }
    }
}

