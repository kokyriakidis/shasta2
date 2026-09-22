#pragma once

// Shasta.
#include "Base.hpp"
#include "ReadId.hpp"

// Standard library.
#include "cstdint.hpp"
#include "string.hpp"
#include "vector.hpp"

namespace shasta2 {
    class Msa1BayesianModel;
}

// The empirical P(observed length m | true length n, base, strand) model
// used by RunLengthEstimator::Bayesian (see msa1.hpp). Built outside
// shasta2, from real reads mapped to a truth reference (see the
// shasta2-homopolymer-model repository's README for the pipeline and the
// exact file format this loads), and loaded here as a small flat table -
// even at generous bounds the whole model is at most a few hundred
// kilobytes, so there is no case for shasta2's mapped-memory machinery.
//
// strand is the read's own OrientedReadId strand (0 = as given by the
// sequencer, 1 = reverse complement) - not the genomic strand a read
// happens to align to, which the model has no way to know and does not
// need to: a homopolymer of base X on strand 1 is definitionally the same
// physical event as one of complement(X) on strand 0, so the model-building
// tool only ever measures strand 0 and fills in strand 1 by complementing
// before writing the file. This class does not know or care about that;
// it just loads whatever four-dimensional table the file contains.
class shasta2::Msa1BayesianModel {
public:

    // The single shared model for a given matrix file, loaded on first use
    // and reused after that. Asserts if called again with a different
    // fileName - in practice a single process run uses exactly one matrix.
    static const Msa1BayesianModel& instance(const string& fileName);

    // log P(m | n, base, strand), Laplace-smoothed (see the .cpp) so this
    // is always finite: an (n, m) combination absent from the training
    // data gets a small but nonzero probability instead of log(0), which
    // would make some candidate n impossible to pick for a reason
    // unrelated to whether it is actually likely. n and m are clamped to
    // this model's bounds (see maxN/maxM), matching the training data's own
    // saturating top bin ("this length or more").
    double logLikelihood(Base base, Strand strand, uint64_t n, uint64_t m) const;

    // log P(n): the training data's own marginal distribution over true
    // length, pooled across base and strand and Laplace-smoothed the same
    // way as logLikelihood, computed once at load time. The simplest
    // defensible default prior; see the .cpp for the rationale and for why
    // a finer granularity (per-base, locus-adaptive) was left for later.
    double logPrior(uint64_t n) const;

    // The largest true length this model has an opinion about (the top,
    // saturating bin). A caller maximizing the posterior over n should
    // search 0..maxN() inclusive, not stop at the largest length actually
    // observed in the column being estimated - unlike every other
    // RunLengthEstimator, Bayesian can legitimately prefer an n no single
    // covering row reported.
    uint64_t maxN() const { return nBins - 1; }
    uint64_t maxM() const { return mBins - 1; }

    // Direct construction from raw counts, bypassing the file format
    // entirely. Used by testMsa1BayesianEstimator() to build small,
    // hand-checkable models without touching the filesystem; production
    // code always goes through instance(fileName) instead. counts is
    // indexed exactly as loadCounts() (see the .cpp) fills it from a file:
    // counts[(((base * 2) + strand) * nBins + n) * mBins + m], row-major,
    // size baseCount(4) * strandCount(2) * nBins * mBins.
    Msa1BayesianModel(uint64_t nBins, uint64_t mBins, const vector<uint64_t>& counts);

private:
    explicit Msa1BayesianModel(const string& fileName);
    void buildFromCounts(const vector<uint64_t>& counts);
    uint64_t cellIndex(uint64_t base, uint64_t strand, uint64_t n, uint64_t m) const;

    uint64_t nBins = 0;
    uint64_t mBins = 0;

    // logLikelihoodTable[cellIndex(base, strand, n, m)].
    vector<double> logLikelihoodTable;

    // logPriorTable[n].
    vector<double> logPriorTable;
};
