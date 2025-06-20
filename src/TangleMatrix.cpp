#include "TangleMatrix.hpp"
#include "Anchor.hpp"
#include "Base.hpp"
using namespace shasta;



TangleMatrix::TangleMatrix(
    const AssemblyGraph& assemblyGraph,
    vector<edge_descriptor> entranceEdges,
    vector<edge_descriptor> exitEdges,
    double aDrift,
    double bDrift)
{

    for(const edge_descriptor e: entranceEdges) {
        const AssemblyGraphEdge& edge = assemblyGraph[e];
        entrances.emplace_back(e, edge.back());
    }
    assemblyGraph.sortEdgeDescriptors(entranceEdges);

    for(const edge_descriptor e: exitEdges) {
        const AssemblyGraphEdge& edge = assemblyGraph[e];
        exits.emplace_back(e, edge.front());
    }
    assemblyGraph.sortEdgeDescriptors(exitEdges);



    // Compute the tangle matrix.
    tangleMatrix.resize(entrances.size(), vector<AnchorPair>(exits.size()));
    for(uint64_t iEntrance=0; iEntrance<entrances.size(); iEntrance++) {
        EntranceOrExit& entrance = entrances[iEntrance];
        for(uint64_t iExit=0; iExit<exits.size(); iExit++) {
            EntranceOrExit& exit = exits[iExit];

            tangleMatrix[iEntrance][iExit] = assemblyGraph.anchors.bridge(
                entrance.step.anchorPair,
                exit.step.anchorPair,
                aDrift, bDrift);
            const uint64_t coverage = tangleMatrix[iEntrance][iExit].orientedReadIds.size();
            entrance.commonCoverage += coverage;
            exit.commonCoverage += coverage;
        }
    }
}



void TangleMatrix::writeHtml(
    const AssemblyGraph& assemblyGraph,
    ostream& html) const
{
    html <<
        "<h2>Tangle matrix</h2>"
        "<p>"
        "<table>"
        "<tr><td rowspan=5 colspan=5>"
        "<th colspan=" << (exits.size() + 1) << " style='background-color:LightCyan'>Exits";

    html << "<tr><th style='background-color:LightCyan'>Index";
    for(uint64_t i=0; i<exits.size(); i++) {
        html << "<td class=centered style='background-color:LightCyan'>" << i;
    }

    html << "<tr><th style='background-color:LightCyan'>Segment";
    for(uint64_t i=0; i<exits.size(); i++) {
        html << "<td class=centered style='background-color:LightCyan'>" << assemblyGraph[exits[i].e].id;
    }

    html << "<tr><th style='background-color:LightCyan'>Coverage";
    for(uint64_t i=0; i<exits.size(); i++) {
        html << "<td class=centered style='background-color:LightCyan'>" << exits[i].coverage();
    }

    html << "<tr><th style='background-color:LightCyan'>Common<br>coverage";
    for(uint64_t i=0; i<exits.size(); i++) {
        html << "<td class=centered style='background-color:LightCyan'>" << exits[i].commonCoverage;
    }

    html <<
        "<tr><th rowspan=" << entrances.size() + 1 <<
        " style='background-color:CornSilk'>Entrances"
        "<th style='background-color:CornSilk'>Index"
        "<th style='background-color:CornSilk'>Segment"
        "<th style='background-color:CornSilk'>Coverage"
        "<th style='background-color:CornSilk'>Common<br>coverage<td>";

    for(uint64_t i=0; i<exits.size(); i++) {
        html << "<td>";
    }

    for(uint64_t iEntrance=0; iEntrance<entrances.size(); iEntrance++) {
        const EntranceOrExit& entrance = entrances[iEntrance];
        html <<
            "<tr>"
            "<td class=centered style='background-color:CornSilk'>" << iEntrance <<
            "<td class=centered style='background-color:CornSilk'>" << assemblyGraph[entrance.e].id <<
            "<td class=centered style='background-color:CornSilk'>" << entrance.coverage() <<
            "<td class=centered style='background-color:CornSilk'>" << entrance.commonCoverage <<
            "<td>";
        for(uint64_t iExit=0; iExit<exits.size(); iExit++) {
            html << "<td class=centered style='background-color:LightPink'>" <<
                tangleMatrix[iEntrance][iExit].orientedReadIds.size();
        }
    }

    html << "</table>";


    if(not hypotheses.empty()) {
        html <<
            "<h3>G test</h3>"
            "<table><tr>"
            "<th>Connectivity<br>matrix"
            "<th>G<br>(dB)";

        for(const auto& hypothesis: hypotheses) {
            html << "<tr><td style='display: flex; align-items: center; justify-content: center;'>";

            html << "<table>";
            for(uint64_t iEntrance=0; iEntrance<entrances.size(); iEntrance++) {
                html << "<tr>";
                for(uint64_t iExit=0; iExit<exits.size(); iExit++) {
                    html << "<td class=centered>" << hypothesis.connectivityMatrix[iEntrance][iExit];
                }
            }
            html << "</table>";

            html <<
                "<td class=centered>" << std::fixed << std::setprecision(1) << hypothesis.G;

        }

        html << "</table>";
    }

}



// Likelihood ratio test of the tangle matrix (G test).
// https://en.wikipedia.org/wiki/G-test
bool TangleMatrix::gTest(double epsilon)
{
    const bool debug = false;

    const uint64_t entranceCount = entrances.size();
    const uint64_t exitCount = exits.size();

    // Limit to a maximum of 16 tangle matrix entries.
    const uint64_t totalTangleMatrixEntryCount = entranceCount * exitCount;
    if(totalTangleMatrixEntryCount > 16) {
        return false;
    }

    // Compute total common coverage.
    uint64_t totalCommonCoverageOnEntrances = 0;
    for(const auto& entrance: entrances) {
        totalCommonCoverageOnEntrances += entrance.commonCoverage;
    }
    uint64_t totalCommonCoverageOnExits = 0;
    for(const auto& exit: exits) {
        totalCommonCoverageOnExits += exit.commonCoverage;
    }
    SHASTA_ASSERT(totalCommonCoverageOnEntrances == totalCommonCoverageOnExits);
    const uint64_t totalCommonCoverage = totalCommonCoverageOnEntrances;

    // Compute what the tangle matrix would be under entirely random assumptions.
    vector< vector<double> > randomTangleMatrix(entranceCount, vector<double>(exitCount, 0.));
    for(uint64_t iEntrance=0; iEntrance<entranceCount; iEntrance++) {
        for(uint64_t iExit=0; iExit<exitCount; iExit++) {
            randomTangleMatrix[iEntrance][iExit] =
                double(entrances[iEntrance].commonCoverage) *
                double(exits[iExit].commonCoverage) /
                double(totalCommonCoverage);
        }
    }

    // Some matrices used in the loop below.
    vector< vector<double> > idealTangleMatrixA(entranceCount, vector<double>(exitCount));
    vector< vector<double> > idealTangleMatrixB(entranceCount, vector<double>(exitCount));
    vector< vector<double> > idealTangleMatrix(entranceCount, vector<double>(exitCount));
    vector< vector<double> > expectedTangleMatrix(entranceCount, vector<double>(exitCount));
    vector< vector<bool> > connectivityMatrix(entranceCount, vector<bool>(exitCount));



    // The connectivity matrix contains true for entrance/exit pairs to be connected.
    // Try all N possible connectivity matrices. Each of them can generate a Hypothesis.
    const uint64_t N = 1ULL << totalTangleMatrixEntryCount;
    hypotheses.clear();
    for(uint64_t connectivityInteger=0; connectivityInteger<N; connectivityInteger++) {

        // Use the bits of connectivityInteger to construct the connectivity matrix.
        uint64_t mask = 1;
        for(uint64_t iEntrance=0; iEntrance<entranceCount; iEntrance++) {
            for(uint64_t iExit=0; iExit<exitCount; iExit++) {
                connectivityMatrix[iEntrance][iExit] = ((connectivityInteger & mask) != 0);
                mask = mask << 1;
            }
        }

        if(debug) {
            cout << "Trying connectivity matrix:" << endl;
            for(uint64_t iEntrance=0; iEntrance<entranceCount; iEntrance++) {
                for(uint64_t iExit=0; iExit<exitCount; iExit++) {
                    cout << connectivityMatrix[iEntrance][iExit] << " ";
                }
                cout << endl;
            }
        }

        // Compute the tangle matrix we would see assuming this connectivity matrix
        // and no errors. This can be done in two ways:
        // - Equally distributing common coverage at each entrance among all the
        //   exits for which the connectivity matrix is true for that entrance (idealTangleMatrixA).
        // - Equally distributing common coverage at each exit among all the
        //   entrances for which the connectivity matrix is true for that entrance (idealTangleMatrixB).
        // We average the two ways to compute the idealTangleMatrix.

        // Compute idealTangleMatrixA.
        for(uint64_t iEntrance=0; iEntrance<entranceCount; iEntrance++) {
            uint64_t nonZeroCount = 0;
            for(uint64_t iExit=0; iExit<exitCount; iExit++) {
                if(connectivityMatrix[iEntrance][iExit]) {
                    ++nonZeroCount;
                }
            }
            if(nonZeroCount == 0) {
                for(uint64_t iExit=0; iExit<exitCount; iExit++) {
                    idealTangleMatrixA[iEntrance][iExit] = 0.;
                }
            } else {
                const double value = double(entrances[iEntrance].commonCoverage) / double(nonZeroCount);
                for(uint64_t iExit=0; iExit<exitCount; iExit++) {
                    if(connectivityMatrix[iEntrance][iExit]) {
                        idealTangleMatrixA[iEntrance][iExit] = value;
                    } else {
                        idealTangleMatrixA[iEntrance][iExit] = 0.;
                    }
                }
            }
        }

        // Compute idealTangleMatrixB.
        for(uint64_t iExit=0; iExit<exitCount; iExit++) {
            uint64_t nonZeroCount = 0;
            for(uint64_t iEntrance=0; iEntrance<entranceCount; iEntrance++) {
                if(connectivityMatrix[iEntrance][iExit]) {
                    ++nonZeroCount;
                }
            }
            if(nonZeroCount == 0) {
                for(uint64_t iEntrance=0; iEntrance<entranceCount; iEntrance++) {
                    idealTangleMatrixB[iEntrance][iExit] = 0.;
                }
            } else {
                const double value = double(exits[iExit].commonCoverage) / double(nonZeroCount);
                for(uint64_t iEntrance=0; iEntrance<entranceCount; iEntrance++) {
                    if(connectivityMatrix[iEntrance][iExit]) {
                        idealTangleMatrixB[iEntrance][iExit] = value;
                    } else {
                        idealTangleMatrixB[iEntrance][iExit] = 0.;
                    }
                }
            }
        }

        // Compute the ideal tangle matrix by averaging idealTangleMatrixA and idealTangleMatrixB.
        for(uint64_t iEntrance=0; iEntrance<entranceCount; iEntrance++) {
            for(uint64_t iExit=0; iExit<exitCount; iExit++) {
                idealTangleMatrix[iEntrance][iExit] =  0.5 *
                    (idealTangleMatrixA[iEntrance][iExit] + idealTangleMatrixB[iEntrance][iExit]);
            }
        }

        // Compute the expected tangle matrix given this connectivity matrix.
        for(uint64_t iEntrance=0; iEntrance<entranceCount; iEntrance++) {
            for(uint64_t iExit=0; iExit<exitCount; iExit++) {
                expectedTangleMatrix[iEntrance][iExit] =
                    (1. - epsilon) * idealTangleMatrix [iEntrance][iExit] +
                    epsilon        * randomTangleMatrix[iEntrance][iExit];
            }
        }

        if(debug) {
            cout << "Expected non-ideal tangle matrix for this connectivity matrix:" << endl;
            for(uint64_t iEntrance=0; iEntrance<entranceCount; iEntrance++) {
                for(uint64_t iExit=0; iExit<exitCount; iExit++) {
                    cout << expectedTangleMatrix[iEntrance][iExit] << " ";
                }
                cout << endl;
            }
        }

        // Do a G-test of the observed tangle matrix against the expected one.
        // https://en.wikipedia.org/wiki/G-test
        double G = 0.;
        for(uint64_t iEntrance=0; iEntrance<entranceCount; iEntrance++) {
            for(uint64_t iExit=0; iExit<exitCount; iExit++) {
                const double actualCoverage = double(tangleMatrix[iEntrance][iExit].size());
                if(actualCoverage > 0.) {
                    const double expectedCoverage = expectedTangleMatrix[iEntrance][iExit];
                    G += actualCoverage * log10(actualCoverage / expectedCoverage);
                }
            }
        }
        G *= 20.;  // Factor of 2 and convert to decibels.

        if(debug) {
            cout << "G " << G << endl;
        }

        // Store this hypothesis.
        hypotheses.emplace_back(Hypothesis(connectivityMatrix, G));
    }

    sort(hypotheses.begin(), hypotheses.end());

    return true;
}



// Return true if there is a single exit for each entrance.
bool TangleMatrix::Hypothesis::isForwardInjective() const
{
    const uint64_t entranceCount = connectivityMatrix.size();
    const uint64_t exitCount = connectivityMatrix.front().size();

    for(uint64_t iEntrance=0; iEntrance<entranceCount; iEntrance++) {
        uint64_t count = 0;
        for(uint64_t iExit=0; iExit<exitCount; iExit++) {
            if(connectivityMatrix[iEntrance][iExit]) {
                ++count;
            }
        }
        if(count != 1) {
            return false;
        }
    }
    return true;
}



// Return true if there is a single entrance for exit entrance.
bool TangleMatrix::Hypothesis::isBackwardInjective() const
{
    const uint64_t entranceCount = connectivityMatrix.size();
    const uint64_t exitCount = connectivityMatrix.front().size();

    for(uint64_t iExit=0; iExit<exitCount; iExit++) {
        uint64_t count = 0;
        for(uint64_t iEntrance=0; iEntrance<entranceCount; iEntrance++) {
            if(connectivityMatrix[iEntrance][iExit]) {
                ++count;
            }
        }
        if(count != 1) {
            return false;
        }
    }
    return true;

}
