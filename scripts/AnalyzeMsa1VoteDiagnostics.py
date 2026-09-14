#!/usr/bin/python3

# Part of the msa1 hard-region evaluation harness. See FindMsa1HardRegions.py
# and EvaluateMsa1AgainstTruth.py, which produce this script's inputs.
#
# For every candidate region where msa1's repair moved the consensus from
# exactly matching HG002 truth to not matching it ("broke"), or from not
# matching to exactly matching ("fixed"), re-runs the repair with
# Assembler.runLocalAssemblyMsa1WithDiagnostics to recover the actual vote
# (Msa1ColumnDiagnostic, see msa1.hpp) at the column where the estimator
# nudged the length, and reports summary statistics for each feature the vote
# carries, comparing the fixed group against the broken group. This is the
# check for whether any feature of the raw vote - not just whichever margin
# or threshold a specific RunLengthEstimator happens to compute from it -
# separates the two groups (i.e. whether a smarter estimator built from the
# same column's vote could do better).

import shasta2

import argparse
import csv
import statistics

parser = argparse.ArgumentParser(description =
    "Compare the actual per-column vote between msa1 regions that got fixed "
    "(wrong to exact) and broke (exact to wrong), against HG002 truth.")
parser.add_argument("stage", type=str,
    help="Assembly stage the candidates were found in (same Data/ directory "
         "used by FindMsa1HardRegions.py).")
parser.add_argument("--candidates", type=str, default="Msa1CandidateRegions.csv")
parser.add_argument("--truth-report", type=str, default="Msa1TruthReport.csv")
arguments = parser.parse_args()

options = shasta2.Options()
assembler = shasta2.Assembler()
assembler.accessAnchors()
assembler.accessJourneys()

candidates = {}
with open(arguments.candidates, newline="") as csvFile:
    for row in csv.DictReader(csvFile):
        candidates[(row["EdgeId"], row["StepIndex"])] = row

fixed = []
broke = []
with open(arguments.truth_report, newline="") as csvFile:
    for row in csv.DictReader(csvFile):
        if not row["TruthStatus"].startswith("ok") or row["DistanceAdaptive"] == "":
            continue
        distanceAdaptive = int(row["DistanceAdaptive"])
        distanceMsa1 = int(row["DistanceMsa1"])
        key = (row["EdgeId"], row["StepIndex"])
        if distanceAdaptive > 0 and distanceMsa1 == 0:
            fixed.append(key)
        elif distanceAdaptive == 0 and distanceMsa1 > 0:
            broke.append(key)
print(len(fixed), "fixed (wrong -> exact),", len(broke), "broke (exact -> wrong)")



def gatherNudgedColumns(keys):
    records = []
    for key in keys:
        candidate = candidates[key]
        anchorIdA = int(candidate["AnchorIdA"])
        anchorIdB = int(candidate["AnchorIdB"])
        orientedReadIds = candidate["OrientedReadIds"].split(";") if candidate["OrientedReadIds"] else []
        success, consensus, diagnostics = assembler.runLocalAssemblyMsa1WithDiagnostics(
            anchorIdA, anchorIdB, orientedReadIds)
        # Keep only the column(s) where the estimator actually nudged the
        # length up by one - that is the decision responsible for the region
        # ending up fixed or broken.
        for (totalWeight, maxObserved, medianLength, cumulativeAtMedian,
                weightAtMedian, weightAtMedianPlusOne, chosenLength) in diagnostics:
            if chosenLength != medianLength + 1:
                continue
            records.append({
                "totalWeight": totalWeight,
                "maxObserved": maxObserved,
                "medianLength": medianLength,
                "margin": cumulativeAtMedian / totalWeight,
                "weightAtMedianShare": weightAtMedian / totalWeight,
                "weightAtMedianPlusOneShare": weightAtMedianPlusOne / totalWeight,
                "rivalRatio": (weightAtMedianPlusOne / weightAtMedian) if weightAtMedian else None,
            })
    return records

fixedRecords = gatherNudgedColumns(fixed)
brokeRecords = gatherNudgedColumns(broke)
print(len(fixedRecords), "nudged columns in fixed regions,",
      len(brokeRecords), "nudged columns in broke regions")



def summarize(records, feature):
    values = sorted(r[feature] for r in records if r[feature] is not None)
    n = len(values)
    if n == 0:
        return None
    return {
        "n": n, "min": values[0], "p25": values[n // 4], "median": values[n // 2],
        "p75": values[3 * n // 4], "max": values[-1], "mean": statistics.mean(values)}

features = [
    "totalWeight", "maxObserved", "medianLength", "margin",
    "weightAtMedianShare", "weightAtMedianPlusOneShare", "rivalRatio"]

print()
print(f"{'feature':<28} {'fixed mean':>12} {'broke mean':>12} {'fixed median':>14} {'broke median':>14}")
for feature in features:
    fixedSummary = summarize(fixedRecords, feature)
    brokeSummary = summarize(brokeRecords, feature)
    if fixedSummary is None or brokeSummary is None:
        continue
    print(f"{feature:<28} {fixedSummary['mean']:>12.3f} {brokeSummary['mean']:>12.3f} "
          f"{fixedSummary['median']:>14.3f} {brokeSummary['median']:>14.3f}")
