// Shasta2.
#include "Msa1BayesianModel.hpp"
#include "SHASTA2_ASSERT.hpp"
using namespace shasta2;

// Standard library.
#include "algorithm.hpp"
#include "fstream.hpp"
#include "stdexcept.hpp"
#include <cmath>



// File format written by shasta2-homopolymer-model/build_matrix.py (see
// that repository's README for the full rationale):
//
//   offset  type      meaning
//   0       char[4]   magic "SHPM"
//   4       uint32    format version (currently 1)
//   8       uint32    baseCount (4 - A,C,G,T)
//   12      uint32    strandCount (2)
//   16      uint32    nBins
//   20      uint32    mBins
//   24      uint64[baseCount][strandCount][nBins][mBins]  raw counts, row-major
//
// Self-describing (dimensions are read from the file, not hardwired here)
// so a future change to the binning in the separate matrix-building
// repository does not require a synchronized hand-edit in this one.
namespace {
    const uint32_t expectedFormatVersion = 1;
    const uint32_t expectedBaseCount = 4;
    const uint32_t expectedStrandCount = 2;
}



uint64_t Msa1BayesianModel::cellIndex(uint64_t base, uint64_t strand, uint64_t n, uint64_t m) const
{
    return (((base * 2) + strand) * nBins + n) * mBins + m;
}



Msa1BayesianModel::Msa1BayesianModel(const string& fileName)
{
    ifstream file(fileName, std::ios::binary);
    if(not file) {
        throw runtime_error("Msa1BayesianModel: could not open " + fileName);
    }

    char magic[4];
    file.read(magic, 4);
    if((not file) or (magic[0] != 'S') or (magic[1] != 'H') or
        (magic[2] != 'P') or (magic[3] != 'M')) {
        throw runtime_error("Msa1BayesianModel: " + fileName +
            " does not have the expected SHPM header.");
    }

    uint32_t formatVersion, baseCount, strandCount, nBins32, mBins32;
    file.read(reinterpret_cast<char*>(&formatVersion), sizeof(formatVersion));
    file.read(reinterpret_cast<char*>(&baseCount), sizeof(baseCount));
    file.read(reinterpret_cast<char*>(&strandCount), sizeof(strandCount));
    file.read(reinterpret_cast<char*>(&nBins32), sizeof(nBins32));
    file.read(reinterpret_cast<char*>(&mBins32), sizeof(mBins32));
    if(not file) {
        throw runtime_error("Msa1BayesianModel: " + fileName + " header is truncated.");
    }
    if(formatVersion != expectedFormatVersion) {
        throw runtime_error("Msa1BayesianModel: " + fileName +
            " has format version " + to_string(formatVersion) +
            ", expected " + to_string(expectedFormatVersion) + ".");
    }
    if((baseCount != expectedBaseCount) or (strandCount != expectedStrandCount)) {
        throw runtime_error("Msa1BayesianModel: " + fileName +
            " has an unexpected base/strand count.");
    }
    nBins = nBins32;
    mBins = mBins32;
    SHASTA2_ASSERT(nBins > 0);
    SHASTA2_ASSERT(mBins > 0);

    const uint64_t cellCount = uint64_t(baseCount) * uint64_t(strandCount) * nBins * mBins;
    vector<uint64_t> counts(cellCount);
    file.read(reinterpret_cast<char*>(counts.data()), std::streamsize(cellCount * sizeof(uint64_t)));
    if(not file) {
        throw runtime_error("Msa1BayesianModel: " + fileName +
            " is truncated - expected " + to_string(cellCount) + " count cells.");
    }

    buildFromCounts(counts);
}



Msa1BayesianModel::Msa1BayesianModel(
    uint64_t nBinsArgument,
    uint64_t mBinsArgument,
    const vector<uint64_t>& counts) :
    nBins(nBinsArgument),
    mBins(mBinsArgument)
{
    SHASTA2_ASSERT(nBins > 0);
    SHASTA2_ASSERT(mBins > 0);
    SHASTA2_ASSERT(counts.size() == 4 * 2 * nBins * mBins);
    buildFromCounts(counts);
}



// Turn raw counts into log-probabilities, with Laplace (add-one) smoothing
// so every cell this model can be asked about has a finite, nonzero
// log-probability - see the header comment on logLikelihood/logPrior for
// why an unsmoothed table (log(0) for an (n, m) never seen in training) is
// not safe to hand to an argmax over n.
void Msa1BayesianModel::buildFromCounts(const vector<uint64_t>& counts)
{
    SHASTA2_ASSERT(counts.size() == 4 * 2 * nBins * mBins);

    logLikelihoodTable.assign(counts.size(), 0.);
    for(uint64_t base=0; base<4; base++) {
        for(uint64_t strand=0; strand<2; strand++) {
            for(uint64_t n=0; n<nBins; n++) {
                double rowTotal = 0.;
                for(uint64_t m=0; m<mBins; m++) {
                    rowTotal += double(counts[cellIndex(base, strand, n, m)]) + 1.;
                }
                for(uint64_t m=0; m<mBins; m++) {
                    const double smoothed = double(counts[cellIndex(base, strand, n, m)]) + 1.;
                    logLikelihoodTable[cellIndex(base, strand, n, m)] = std::log(smoothed / rowTotal);
                }
            }
        }
    }

    // The prior: the marginal over n, pooled across base and strand, from
    // the same raw counts, with the same Laplace smoothing.
    vector<double> priorCount(nBins, 0.);
    for(uint64_t base=0; base<4; base++) {
        for(uint64_t strand=0; strand<2; strand++) {
            for(uint64_t n=0; n<nBins; n++) {
                for(uint64_t m=0; m<mBins; m++) {
                    priorCount[n] += double(counts[cellIndex(base, strand, n, m)]);
                }
            }
        }
    }
    double priorTotal = 0.;
    for(uint64_t n=0; n<nBins; n++) {
        priorCount[n] += 1.;
        priorTotal += priorCount[n];
    }
    logPriorTable.assign(nBins, 0.);
    for(uint64_t n=0; n<nBins; n++) {
        logPriorTable[n] = std::log(priorCount[n] / priorTotal);
    }
}



double Msa1BayesianModel::logLikelihood(Base base, Strand strand, uint64_t n, uint64_t m) const
{
    const uint64_t nClamped = min(n, maxN());
    const uint64_t mClamped = min(m, maxM());
    return logLikelihoodTable[cellIndex(base.value, strand, nClamped, mClamped)];
}



double Msa1BayesianModel::logPrior(uint64_t n) const
{
    return logPriorTable[min(n, maxN())];
}



const Msa1BayesianModel& Msa1BayesianModel::instance(const string& fileName)
{
    // Thread-safe by the language's own guarantee on static local
    // initialization ("magic statics") - no explicit locking needed here.
    // A process is only ever expected to use one Bayesian matrix, so
    // requesting a different fileName after the first call is treated as a
    // programming error rather than silently ignored or reloaded.
    static const string loadedFileName = fileName;
    static const Msa1BayesianModel model(fileName);

    if(fileName != loadedFileName) {
        throw runtime_error("Msa1BayesianModel::instance called with " + fileName +
            " after already having loaded " + loadedFileName +
            " - a process can only use one Bayesian matrix.");
    }
    return model;
}
