#!/usr/bin/python3

# Part of the msa1 hard-region evaluation harness.
#
# Walks every AssemblyGraph step of the given assembly stage and re-runs
# LocalAssembly7 twice per step, both with Method::Adaptive, once with
# Options::useMsa1 false and once true. Most steps produce an identical
# consensus either way - msa1's own internal trigger and impure-column
# safety net didn't fire or didn't change anything. This writes a CSV row
# for every step where the two consensuses differ, i.e. every step where
# the repair actually did something. That CSV is the input to
# EvaluateMsa1AgainstTruth.py, which checks whether the repair moved the
# consensus closer to or further from HG002 truth.

import shasta2

import argparse
import csv

parser = argparse.ArgumentParser(description =
    "Find AssemblyGraph steps where msa1's homopolymer repair changes the consensus.")
parser.add_argument("stage", type=str, help="Assembly stage to evaluate.")
parser.add_argument("--output", type=str, default="Msa1CandidateRegions.csv",
    help="Output csv file name.")
parser.add_argument("--max-steps", type=int, default=None,
    help="Only examine the first this many steps (for a quick smoke test).")
arguments = parser.parse_args()

options = shasta2.Options()

assembler = shasta2.Assembler()
assembler.accessAnchors()
assembler.accessJourneys()

assemblyGraph = assembler.getAssemblyGraph(arguments.stage, options)

steps = assemblyGraph.getAssemblyGraphSteps()
if arguments.max_steps is not None:
    steps = steps[:arguments.max_steps]
print(len(steps), "AssemblyGraph steps to examine.")

candidateCount = 0
failureCount = 0
with open(arguments.output, "w", newline="") as csvFile:
    writer = csv.writer(csvFile)
    writer.writerow([
        "EdgeId", "StepIndex", "AnchorIdA", "AnchorIdB",
        "OrientedReadIds", "ConsensusNoRepair", "ConsensusWithRepair"])

    for i, (edgeId, stepIndex, anchorIdA, anchorIdB, orientedReadIdStrings) in enumerate(steps):
        if i % 10000 == 0 and i > 0:
            print(i, "steps examined,", candidateCount, "candidates found so far.")

        successNoRepair, consensusNoRepair, successWithRepair, consensusWithRepair = \
            assembler.runLocalAssemblyWithAndWithoutMsa1Repair(anchorIdA, anchorIdB, orientedReadIdStrings)

        if not (successNoRepair and successWithRepair):
            failureCount += 1
            continue

        if consensusNoRepair == consensusWithRepair:
            continue

        candidateCount += 1
        writer.writerow([
            edgeId, stepIndex, anchorIdA, anchorIdB,
            ";".join(orientedReadIdStrings),
            consensusNoRepair, consensusWithRepair])

print(len(steps), "steps examined.")
print(candidateCount, "candidate hard regions found (msa1 repair changed the consensus).")
print(failureCount, "steps failed to assemble (skipped) with or without the repair.")
