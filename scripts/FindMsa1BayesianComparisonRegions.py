#!/usr/bin/python3

# Compares RunLengthEstimator::MedianMarginGated (the production default)
# against RunLengthEstimator::Bayesian directly, on the same real assembly.
#
# Scans every AssemblyGraph step in parallel (see
# AssemblyGraph::findMsa1BayesianComparisonRegions), running LocalAssembly7
# twice per step - once at each estimator - and writes a CSV row for every
# step where EITHER estimator's repair changed the consensus relative to no
# repair. That CSV is the input to EvaluateMsa1BayesianAgainstTruth.py, which
# checks which estimator (if either) moved the consensus closer to HG002
# truth.

import shasta2

import argparse
import csv
import time

parser = argparse.ArgumentParser(description =
    "Find AssemblyGraph steps where msa1's MedianMarginGated and/or Bayesian "
    "repair changes the consensus, for comparing the two directly.")
parser.add_argument("stage", type=str, help="Assembly stage to evaluate.")
parser.add_argument("bayesian_matrix", type=str,
    help="Path to the Bayesian estimator's P(m|n,base,strand) matrix file.")
parser.add_argument("--output", type=str, default="Msa1BayesianComparisonRegions.csv",
    help="Output csv file name.")
arguments = parser.parse_args()

options = shasta2.Options()

assembler = shasta2.Assembler()
assembler.accessAnchors()
assembler.accessJourneys()

assemblyGraph = assembler.getAssemblyGraph(arguments.stage, options)

t0 = time.time()
candidates = assemblyGraph.findMsa1BayesianComparisonRegions(arguments.bayesian_matrix)
t1 = time.time()
print(f"Scan completed in {t1 - t0:.1f} seconds.")
print(len(candidates), "candidate regions found (either estimator's repair changed the consensus).")

with open(arguments.output, "w", newline="") as csvFile:
    writer = csv.writer(csvFile)
    writer.writerow([
        "EdgeId", "StepIndex", "AnchorIdA", "AnchorIdB",
        "OrientedReadIds", "ConsensusNoRepair", "ConsensusMedian", "ConsensusBayesian"])

    for edgeId, stepIndex, anchorIdA, anchorIdB, orientedReadIdStrings, \
            consensusNoRepair, consensusMedian, consensusBayesian in candidates:
        writer.writerow([
            edgeId, stepIndex, anchorIdA, anchorIdB,
            ";".join(orientedReadIdStrings),
            consensusNoRepair, consensusMedian, consensusBayesian])

print("Written to", arguments.output)
