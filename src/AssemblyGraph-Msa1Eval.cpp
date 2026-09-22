// Shasta2.
#include "AssemblyGraph.hpp"
#include "deduplicate.hpp"
#include "LocalAssembly7.hpp"
#include "msa1.hpp"
#include "Options.hpp"
#include "performanceLog.hpp"
#include "timestamp.hpp"
using namespace shasta2;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>

// Standard library.
#include "algorithm.hpp"
#include "iterator.hpp"
#include <numeric>
#include "utility.hpp"



namespace shasta2 {

    // The oriented reads a step's LocalAssembly7 run is given: the ones on
    // the step's own AnchorPair, plus the ones borrowed from the previous/next
    // step (which is why an oriented read here is not guaranteed to actually
    // be on this step's anchorIdA/anchorIdB - see Assembler::anchorContainsOrientedRead).
    // Shared by getAssemblyGraphSteps and findMsa1CandidateRegions so the two
    // cannot drift apart on what a "step" means.
    static vector<OrientedReadId> msa1StepOrientedReadIds(
        const AssemblyGraphEdge& edge,
        uint64_t i)
    {
        vector<OrientedReadId> additionalOrientedReadIds;
        if(i > 0) {
            std::ranges::copy(edge[i - 1].anchorPair.orientedReadIds, back_inserter(additionalOrientedReadIds));
        }
        if(i < edge.size() - 1) {
            std::ranges::copy(edge[i + 1].anchorPair.orientedReadIds, back_inserter(additionalOrientedReadIds));
        }
        deduplicate(additionalOrientedReadIds);

        vector<OrientedReadId> orientedReadIds = additionalOrientedReadIds;
        std::ranges::copy(edge[i].anchorPair.orientedReadIds, back_inserter(orientedReadIds));
        deduplicate(orientedReadIds);
        return orientedReadIds;
    }
}



// For the msa1 hard-region evaluation harness. Return, for every step of every
// edge, exactly the (anchorIdA, anchorIdB, orientedReadIds) that assembleStep
// would pass to LocalAssembly7 - including the additional oriented reads
// borrowed from the previous/next step. See assembleStep for the logic mirrored here.
vector< tuple<uint64_t, uint64_t, AnchorId, AnchorId, vector<string> > >
    AssemblyGraph::getAssemblyGraphSteps() const
{
    const AssemblyGraph& assemblyGraph = *this;

    vector< tuple<uint64_t, uint64_t, AnchorId, AnchorId, vector<string> > > steps;
    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        const AssemblyGraphEdge& edge = assemblyGraph[e];

        for(uint64_t i=0; i<edge.size(); i++) {
            const vector<OrientedReadId> orientedReadIds = msa1StepOrientedReadIds(edge, i);

            vector<string> orientedReadIdStrings;
            orientedReadIdStrings.reserve(orientedReadIds.size());
            for(const OrientedReadId orientedReadId: orientedReadIds) {
                orientedReadIdStrings.push_back(orientedReadId.getString());
            }

            const AnchorPair& anchorPair = edge[i].anchorPair;
            steps.push_back(make_tuple(
                uint64_t(edge.id), i,
                anchorPair.anchorIdA, anchorPair.anchorIdB,
                orientedReadIdStrings));
        }
    }

    return steps;
}



// See AssemblyGraph.hpp for comments.
vector< tuple<uint64_t, uint64_t, AnchorId, AnchorId, vector<string>, string, string> >
    AssemblyGraph::findMsa1CandidateRegions()
{
    AssemblyGraph& assemblyGraph = *this;

    msa1StepsToScan.clear();
    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        const AssemblyGraphEdge& edge = assemblyGraph[e];
        for(uint64_t i=0; i<edge.size(); i++) {
            msa1StepsToScan.push_back(make_pair(e, i));
        }
    }
    performanceLog << timestamp << "findMsa1CandidateRegions: scanning " <<
        msa1StepsToScan.size() << " steps." << endl;

    const uint64_t threadCount = options.actualThreadCount();
    msa1CandidatesByThread.clear();
    msa1CandidatesByThread.resize(threadCount);

    const uint64_t batchCount = 1;
    setupLoadBalancing(msa1StepsToScan.size(), batchCount);
    runThreads(&AssemblyGraph::findMsa1CandidateRegionsThreadFunction, threadCount);

    // Flatten the per-thread results. Each thread only ever wrote to its own
    // slot, so this is the only part of the whole scan that touches all of
    // them at once.
    uint64_t candidateCount = 0;
    for(const auto& v: msa1CandidatesByThread) {
        candidateCount += v.size();
    }
    vector< tuple<uint64_t, uint64_t, AnchorId, AnchorId, vector<string>, string, string> > candidates;
    candidates.reserve(candidateCount);
    for(auto& v: msa1CandidatesByThread) {
        for(auto& row: v) {
            candidates.push_back(std::move(row));
        }
    }
    msa1CandidatesByThread.clear();
    msa1CandidatesByThread.shrink_to_fit();
    msa1StepsToScan.clear();
    msa1StepsToScan.shrink_to_fit();

    performanceLog << timestamp << "findMsa1CandidateRegions: found " <<
        candidates.size() << " candidate regions." << endl;

    return candidates;
}



void AssemblyGraph::findMsa1CandidateRegionsThreadFunction(uint64_t threadId)
{
    const AssemblyGraph& assemblyGraph = *this;
    ostream html(0);

    vector< tuple<uint64_t, uint64_t, AnchorId, AnchorId, vector<string>, string, string> >&
        candidates = msa1CandidatesByThread[threadId];

    LocalAssembly7::Options localAssembly7Options;
    localAssembly7Options.useMsa1 = true;

    uint64_t begin, end;
    while(getNextBatch(begin, end)) {
        for(uint64_t j=begin; j!=end; j++) {

            if((j % 100000) == 0) {
                std::lock_guard<std::mutex> lock(mutex);
                performanceLog << timestamp << "findMsa1CandidateRegions: step " << j <<
                    " of " << msa1StepsToScan.size() << endl;
            }

            const auto& p = msa1StepsToScan[j];
            const edge_descriptor e = p.first;
            const uint64_t i = p.second;
            const AssemblyGraphEdge& edge = assemblyGraph[e];
            const AnchorPair& anchorPair = edge[i].anchorPair;
            const vector<OrientedReadId> orientedReadIds = msa1StepOrientedReadIds(edge, i);

            // A step that fails to assemble is skipped, not fatal: unlike
            // assembleThreadFunction, this is an evaluation scan over the
            // whole assembly, not the assembly itself, and a handful of
            // difficult steps should not abort a multi-hour run.
            bool success = false;
            vector<shasta2::Base> consensusNoRepair, consensusWithRepair;
            try {
                const LocalAssembly7 localAssembly(
                    localAssembly7Options, anchors,
                    anchorPair.anchorIdA, anchorPair.anchorIdB,
                    html, orientedReadIds);
                success = localAssembly.success;
                consensusNoRepair = localAssembly.sequenceBeforeRepair;
                consensusWithRepair = localAssembly.sequence;
            } catch(std::exception&) {
                success = false;
            }
            if(not success) {
                continue;
            }
            if(consensusNoRepair == consensusWithRepair) {
                continue;
            }

            vector<string> orientedReadIdStrings;
            orientedReadIdStrings.reserve(orientedReadIds.size());
            for(const OrientedReadId orientedReadId: orientedReadIds) {
                orientedReadIdStrings.push_back(orientedReadId.getString());
            }

            candidates.push_back(make_tuple(
                uint64_t(edge.id), i, anchorPair.anchorIdA, anchorPair.anchorIdB,
                orientedReadIdStrings, toString(consensusNoRepair), toString(consensusWithRepair)));
        }
    }
}



// See AssemblyGraph.hpp for comments.
vector< tuple<uint64_t, uint64_t, AnchorId, AnchorId, vector<string>, string, string, string> >
    AssemblyGraph::findMsa1BayesianComparisonRegions(const string& bayesianMatrixName)
{
    AssemblyGraph& assemblyGraph = *this;

    msa1BayesianMatrixNameForScan = bayesianMatrixName;
    msa1BayesianStepsToScan.clear();
    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        const AssemblyGraphEdge& edge = assemblyGraph[e];
        for(uint64_t i=0; i<edge.size(); i++) {
            msa1BayesianStepsToScan.push_back(make_pair(e, i));
        }
    }
    performanceLog << timestamp << "findMsa1BayesianComparisonRegions: scanning " <<
        msa1BayesianStepsToScan.size() << " steps." << endl;

    const uint64_t threadCount = options.actualThreadCount();
    msa1BayesianCandidatesByThread.clear();
    msa1BayesianCandidatesByThread.resize(threadCount);

    const uint64_t batchCount = 1;
    setupLoadBalancing(msa1BayesianStepsToScan.size(), batchCount);
    runThreads(&AssemblyGraph::findMsa1BayesianComparisonThreadFunction, threadCount);

    uint64_t candidateCount = 0;
    for(const auto& v: msa1BayesianCandidatesByThread) {
        candidateCount += v.size();
    }
    vector< tuple<uint64_t, uint64_t, AnchorId, AnchorId, vector<string>, string, string, string> > candidates;
    candidates.reserve(candidateCount);
    for(auto& v: msa1BayesianCandidatesByThread) {
        for(auto& row: v) {
            candidates.push_back(std::move(row));
        }
    }
    msa1BayesianCandidatesByThread.clear();
    msa1BayesianCandidatesByThread.shrink_to_fit();
    msa1BayesianStepsToScan.clear();
    msa1BayesianStepsToScan.shrink_to_fit();

    performanceLog << timestamp << "findMsa1BayesianComparisonRegions: found " <<
        candidates.size() << " candidate regions." << endl;

    return candidates;
}



void AssemblyGraph::findMsa1BayesianComparisonThreadFunction(uint64_t threadId)
{
    const AssemblyGraph& assemblyGraph = *this;
    ostream html(0);

    vector< tuple<uint64_t, uint64_t, AnchorId, AnchorId, vector<string>, string, string, string> >&
        candidates = msa1BayesianCandidatesByThread[threadId];

    LocalAssembly7::Options medianOptions;
    medianOptions.useMsa1 = true;
    medianOptions.estimator = RunLengthEstimator::MedianMarginGated;

    LocalAssembly7::Options bayesianOptions;
    bayesianOptions.useMsa1 = true;
    bayesianOptions.estimator = RunLengthEstimator::Bayesian;
    bayesianOptions.msa1BayesianMatrixName = msa1BayesianMatrixNameForScan;

    uint64_t begin, end;
    while(getNextBatch(begin, end)) {
        for(uint64_t j=begin; j!=end; j++) {

            if((j % 100000) == 0) {
                std::lock_guard<std::mutex> lock(mutex);
                performanceLog << timestamp << "findMsa1BayesianComparisonRegions: step " << j <<
                    " of " << msa1BayesianStepsToScan.size() << endl;
            }

            const auto& p = msa1BayesianStepsToScan[j];
            const edge_descriptor e = p.first;
            const uint64_t i = p.second;
            const AssemblyGraphEdge& edge = assemblyGraph[e];
            const AnchorPair& anchorPair = edge[i].anchorPair;
            const vector<OrientedReadId> orientedReadIds = msa1StepOrientedReadIds(edge, i);

            bool medianSuccess = false, bayesianSuccess = false;
            vector<shasta2::Base> consensusNoRepair, consensusMedian, consensusBayesian;
            try {
                const LocalAssembly7 medianAssembly(
                    medianOptions, anchors,
                    anchorPair.anchorIdA, anchorPair.anchorIdB,
                    html, orientedReadIds);
                medianSuccess = medianAssembly.success;
                consensusNoRepair = medianAssembly.sequenceBeforeRepair;
                consensusMedian = medianAssembly.sequence;
            } catch(std::exception&) {
                medianSuccess = false;
            }
            if(not medianSuccess) {
                continue;
            }

            try {
                const LocalAssembly7 bayesianAssembly(
                    bayesianOptions, anchors,
                    anchorPair.anchorIdA, anchorPair.anchorIdB,
                    html, orientedReadIds);
                bayesianSuccess = bayesianAssembly.success;
                consensusBayesian = bayesianAssembly.sequence;
            } catch(std::exception&) {
                bayesianSuccess = false;
            }
            if(not bayesianSuccess) {
                continue;
            }

            // A candidate step is one where EITHER estimator's repair
            // changed anything relative to no repair - see the comment on
            // this function's declaration.
            if((consensusNoRepair == consensusMedian) and (consensusNoRepair == consensusBayesian)) {
                continue;
            }

            vector<string> orientedReadIdStrings;
            orientedReadIdStrings.reserve(orientedReadIds.size());
            for(const OrientedReadId orientedReadId: orientedReadIds) {
                orientedReadIdStrings.push_back(orientedReadId.getString());
            }

            candidates.push_back(make_tuple(
                uint64_t(edge.id), i, anchorPair.anchorIdA, anchorPair.anchorIdB,
                orientedReadIdStrings, toString(consensusNoRepair),
                toString(consensusMedian), toString(consensusBayesian)));
        }
    }
}



// See AssemblyGraph.hpp for comments.
int64_t shasta2::editDistance(const string& a, const string& b, uint64_t cap)
{
    if(uint64_t(a.size()) * uint64_t(b.size()) > cap) {
        return -1;
    }

    vector<uint32_t> previousRow(b.size() + 1);
    std::iota(previousRow.begin(), previousRow.end(), 0);
    vector<uint32_t> currentRow(b.size() + 1);

    for(uint64_t i=1; i<=a.size(); i++) {
        currentRow[0] = uint32_t(i);
        const char ca = a[i - 1];
        for(uint64_t j=1; j<=b.size(); j++) {
            currentRow[j] = std::min({
                previousRow[j] + 1,
                currentRow[j - 1] + 1,
                previousRow[j - 1] + uint32_t(ca != b[j - 1])});
        }
        std::swap(previousRow, currentRow);
    }

    return int64_t(previousRow[b.size()]);
}
