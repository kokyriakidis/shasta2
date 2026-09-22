#!/usr/bin/python3

# Part of the msa1 hard-region evaluation harness.
#
# Scans every AssemblyGraph step of the given assembly stage, in parallel
# (see AssemblyGraph::findMsa1CandidateRegions), running LocalAssembly7 once
# per step with Method::Adaptive and Options::useMsa1, and comparing the
# consensus before and after the repair. Most steps produce an identical
# consensus either way - msa1's own internal trigger and impure-column
# safety net didn't fire or didn't change anything. This writes a CSV row
# for every step where the two consensuses differ, i.e. every step where
# the repair actually did something. That CSV is the input to
# EvaluateMsa1AgainstTruth.py, which checks whether the repair moved the
# consensus closer to or further from HG002 truth.

import shasta2

import argparse
import csv
import time

parser = argparse.ArgumentParser(description =
    "Find AssemblyGraph steps where msa1's homopolymer repair changes the consensus.")
parser.add_argument("stage", type=str, help="Assembly stage to evaluate.")
parser.add_argument("--output", type=str, default="Msa1CandidateRegions.csv",
    help="Output csv file name.")
arguments = parser.parse_args()

options = shasta2.Options()

assembler = shasta2.Assembler()
assembler.accessAnchors()
assembler.accessJourneys()

assemblyGraph = assembler.getAssemblyGraph(arguments.stage, options)

t0 = time.time()
candidates = assemblyGraph.findMsa1CandidateRegions()
t1 = time.time()
print(f"Scan completed in {t1 - t0:.1f} seconds.")
print(len(candidates), "candidate hard regions found (msa1 repair changed the consensus).")

with open(arguments.output, "w", newline="") as csvFile:
    writer = csv.writer(csvFile)
    writer.writerow([
        "EdgeId", "StepIndex", "AnchorIdA", "AnchorIdB",
        "OrientedReadIds", "ConsensusNoRepair", "ConsensusWithRepair"])

    for edgeId, stepIndex, anchorIdA, anchorIdB, orientedReadIdStrings, \
            consensusNoRepair, consensusWithRepair in candidates:
        writer.writerow([
            edgeId, stepIndex, anchorIdA, anchorIdB,
            ";".join(orientedReadIdStrings),
            consensusNoRepair, consensusWithRepair])

print("Written to", arguments.output)
