#!/usr/bin/python3

# Compares RunLengthEstimator::MedianMarginGated (the production default)
# against RunLengthEstimator::Bayesian against HG002 truth, for the
# candidate regions found by FindMsa1BayesianComparisonRegions.py.
#
# Ground truth establishment is identical to EvaluateMsa1AgainstTruth.py
# (windowed read mapping, MAPQ>=10 gate, CIGAR-based coordinate lifting,
# majority vote across core reads) - see that script's comments for the
# full rationale and the measurements behind those choices. The only
# difference is comparing THREE consensus columns (no repair, Median,
# Bayesian) against the same truth substring instead of two.

import shasta2

import argparse
import csv
import os
import subprocess
import tempfile
from collections import defaultdict



def reverseComplement(sequence):
    complement = str.maketrans("ACGT", "TGCA")
    return sequence.translate(complement)[::-1]



def editDistance(a, b, cap=4_000_000):
    distance = shasta2.editDistance(a, b, cap)
    return None if distance < 0 else distance



QUERY_CONSUMING = set("MIS=X")
REFERENCE_CONSUMING = set("MDN=X")

def parseCigar(cigar):
    ops = []
    length = 0
    for c in cigar:
        if c.isdigit():
            length = length * 10 + int(c)
        else:
            ops.append((length, c))
            length = 0
    return ops



def queryToReferencePosition(cigar, alignmentStart, queryPosition):
    q = 0
    r = alignmentStart
    for length, op in parseCigar(cigar):
        queryConsumes = op in QUERY_CONSUMING
        referenceConsumes = op in REFERENCE_CONSUMING
        if q <= queryPosition < q + length:
            return (r + (queryPosition - q)) if referenceConsumes else None
        if queryConsumes:
            q += length
        if referenceConsumes:
            r += length
    return None



parser = argparse.ArgumentParser(description =
    "Compare msa1's MedianMarginGated- and Bayesian-repaired consensus against "
    "HG002 v1.1 truth for the candidate regions found by "
    "FindMsa1BayesianComparisonRegions.py.")
parser.add_argument("stage", type=str,
    help="Assembly stage the candidates were found in (must be reachable with the "
         "same Data/ directory used by FindMsa1BayesianComparisonRegions.py).")
parser.add_argument("--candidates", type=str, default="Msa1BayesianComparisonRegions.csv")
parser.add_argument("--truth-fasta", type=str,
    default=os.path.expanduser("~/Downloads/hg002v1.1.fasta"))
parser.add_argument("--minimap2", type=str,
    default=os.path.expanduser("~/Downloads/minimap2-2.30_x64-linux/minimap2"))
parser.add_argument("--samtools", type=str, default="samtools")
parser.add_argument("--truth-index", type=str, default=None)
parser.add_argument("--output", type=str, default="Msa1BayesianTruthReport.csv")
parser.add_argument("--work-dir", type=str, default=None)
parser.add_argument("--window-pad", type=int, default=5000)
parser.add_argument("--min-mapq", type=int, default=10)
parser.add_argument("--threads", type=int, default=os.cpu_count() or 4)
arguments = parser.parse_args()

options = shasta2.Options()
assembler = shasta2.Assembler()
assembler.accessAnchors()
assembler.accessJourneys()

truthIndex = arguments.truth_index or (arguments.truth_fasta + ".map-ont.mmi")
if not os.path.exists(truthIndex):
    print("Building minimap2 index", truthIndex, "(one-time cost for this truth fasta)...")
    subprocess.run(
        [arguments.minimap2, "-x", "map-ont", "-d", truthIndex, arguments.truth_fasta],
        check=True)



rows = []
with open(arguments.candidates, newline="") as csvFile:
    for row in csv.DictReader(csvFile):
        row["OrientedReadIds"] = row["OrientedReadIds"].split(";") if row["OrientedReadIds"] else []
        rows.append(row)
print(len(rows), "candidate regions loaded from", arguments.candidates)



for row in rows:
    anchorIdA = int(row["AnchorIdA"])
    anchorIdB = int(row["AnchorIdB"])
    row["CoreReads"] = [
        orientedReadId for orientedReadId in row["OrientedReadIds"]
        if assembler.anchorContainsOrientedRead(anchorIdA, orientedReadId)
        and assembler.anchorContainsOrientedRead(anchorIdB, orientedReadId)]

allReads = sorted({orientedReadId for row in rows for orientedReadId in row["CoreReads"]})
readSequences = {
    orientedReadId: assembler.getOrientedReadSequenceString(orientedReadId)
    for orientedReadId in allReads}
print(len(allReads), "distinct oriented reads among the candidate regions.")



workDir = arguments.work_dir or tempfile.mkdtemp(prefix="msa1BayesianEval_")
os.makedirs(workDir, exist_ok=True)

windowInfo = {}
readsFastaName = os.path.join(workDir, "reads.fasta")
with open(readsFastaName, "w") as fasta:
    for rowIndex, row in enumerate(rows):
        anchorIdA = int(row["AnchorIdA"])
        anchorIdB = int(row["AnchorIdB"])
        for orientedReadId in row["CoreReads"]:
            sequence = readSequences[orientedReadId]
            positionA = assembler.getAnchorPositionInOrientedRead(anchorIdA, orientedReadId)
            positionB = assembler.getAnchorPositionInOrientedRead(anchorIdB, orientedReadId)
            windowBegin = max(0, min(positionA, positionB) - arguments.window_pad)
            windowEnd = min(len(sequence), max(positionA, positionB) + arguments.window_pad)
            windowSequence = sequence[windowBegin:windowEnd]
            windowInfo[(rowIndex, orientedReadId)] = (
                positionA - windowBegin, positionB - windowBegin, len(windowSequence))
            fasta.write(f">{rowIndex}_{orientedReadId}\n{windowSequence}\n")
print(len(windowInfo), "region/read windows to map to truth.")

samFileName = os.path.join(workDir, "reads.sam")
with open(samFileName, "w") as samFile:
    subprocess.run(
        [arguments.minimap2, "-a", "--eqx", "-x", "map-ont", "-t", str(arguments.threads),
            truthIndex, readsFastaName],
        stdout=samFile, stderr=subprocess.DEVNULL, check=True)

alignments = {}
with open(samFileName) as samFile:
    for line in samFile:
        if line.startswith("@"):
            continue
        fields = line.rstrip("\n").split("\t")
        queryName = fields[0]
        flag = int(fields[1])
        if flag & 0x904:
            continue
        referenceName = fields[2]
        referenceStart = int(fields[3]) - 1
        mapq = int(fields[4])
        cigar = fields[5]
        alignments[queryName] = (referenceName, referenceStart, cigar, bool(flag & 16), mapq)
trustedAlignments = sum(1 for a in alignments.values() if a[4] >= arguments.min_mapq)
print(len(alignments), "of", len(windowInfo), "windows have a primary mapping to truth,",
    trustedAlignments, f"at MAPQ>={arguments.min_mapq} or better.")



def fetchTruthSubstring(referenceName, begin, end):
    region = f"{referenceName}:{begin + 1}-{end}"
    result = subprocess.run(
        [arguments.samtools, "faidx", arguments.truth_fasta, region],
        capture_output=True, text=True, check=True)
    return "".join(result.stdout.splitlines()[1:])



# For each region, establish truth once, then score all three consensuses
# against it. verdict compares Bayesian against Median directly (not
# against no-repair) - that is the actual question this experiment asks.
counts = defaultdict(int)

with open(arguments.output, "w", newline="") as outputFile:
    writer = csv.writer(outputFile)
    writer.writerow([
        "EdgeId", "StepIndex", "AnchorIdA", "AnchorIdB",
        "TruthStatus", "TruthLength",
        "DistanceNoRepair", "DistanceMedian", "DistanceBayesian",
        "MedianVerdict", "BayesianVerdict", "BayesianVsMedian",
        "Truth", "ConsensusNoRepair", "ConsensusMedian", "ConsensusBayesian"])

    for rowIndex, row in enumerate(rows):
        anchorIdA = int(row["AnchorIdA"])
        anchorIdB = int(row["AnchorIdB"])

        intervals = []
        for orientedReadId in row["CoreReads"]:
            alignment = alignments.get(f"{rowIndex}_{orientedReadId}")
            if alignment is None:
                continue
            referenceName, referenceStart, cigar, isReverse, mapq = alignment
            if mapq < arguments.min_mapq:
                continue

            localPositionA, localPositionB, windowLength = windowInfo[(rowIndex, orientedReadId)]
            queryPositionA = (windowLength - 1 - localPositionA) if isReverse else localPositionA
            queryPositionB = (windowLength - 1 - localPositionB) if isReverse else localPositionB

            referencePositionA = queryToReferencePosition(cigar, referenceStart, queryPositionA)
            referencePositionB = queryToReferencePosition(cigar, referenceStart, queryPositionB)
            if referencePositionA is None or referencePositionB is None:
                continue

            if isReverse:
                begin = min(referencePositionA, referencePositionB) + 1
                end = max(referencePositionA, referencePositionB) + 1
            else:
                begin = min(referencePositionA, referencePositionB)
                end = max(referencePositionA, referencePositionB)
            if begin >= end:
                continue
            intervals.append((referenceName, begin, end, isReverse))

        consensusNoRepair = row["ConsensusNoRepair"]
        consensusMedian = row["ConsensusMedian"]
        consensusBayesian = row["ConsensusBayesian"]

        if not intervals:
            writer.writerow([
                row["EdgeId"], row["StepIndex"], anchorIdA, anchorIdB,
                "no_reads_mapped", "", "", "", "", "skipped", "skipped", "skipped",
                "", consensusNoRepair, consensusMedian, consensusBayesian])
            counts["skipped"] += 1
            continue

        groups = defaultdict(list)
        for referenceName, begin, end, isReverse in intervals:
            groups[(referenceName, isReverse)].append((begin, end))
        (referenceName, isReverse), bestList = max(groups.items(), key=lambda kv: len(kv[1]))

        if len(bestList) < max(1, len(intervals) // 2 + 1):
            writer.writerow([
                row["EdgeId"], row["StepIndex"], anchorIdA, anchorIdB,
                f"reads_disagree:{len(bestList)}/{len(intervals)}",
                "", "", "", "", "skipped", "skipped", "skipped",
                "", consensusNoRepair, consensusMedian, consensusBayesian])
            counts["skipped"] += 1
            continue

        begins = sorted(b for b, e in bestList)
        ends = sorted(e for b, e in bestList)
        begin = begins[len(begins) // 2]
        end = ends[len(ends) // 2]

        truthSubstring = fetchTruthSubstring(referenceName, begin, end)
        if isReverse:
            truthSubstring = reverseComplement(truthSubstring)

        distanceNoRepair = editDistance(consensusNoRepair, truthSubstring)
        distanceMedian = editDistance(consensusMedian, truthSubstring)
        distanceBayesian = editDistance(consensusBayesian, truthSubstring)

        if None in (distanceNoRepair, distanceMedian, distanceBayesian):
            writer.writerow([
                row["EdgeId"], row["StepIndex"], anchorIdA, anchorIdB,
                f"ok:{len(bestList)}/{len(intervals)}:{referenceName}", len(truthSubstring),
                distanceNoRepair, distanceMedian, distanceBayesian,
                "skipped_too_long", "skipped_too_long", "skipped_too_long",
                truthSubstring, consensusNoRepair, consensusMedian, consensusBayesian])
            counts["skipped_too_long"] += 1
            continue

        def verdict(distanceRepaired):
            if distanceRepaired < distanceNoRepair:
                return "helped"
            if distanceRepaired > distanceNoRepair:
                return "hurt"
            return "neutral"
        medianVerdict = verdict(distanceMedian)
        bayesianVerdict = verdict(distanceBayesian)

        if distanceBayesian < distanceMedian:
            bayesianVsMedian = "bayesian_better"
        elif distanceBayesian > distanceMedian:
            bayesianVsMedian = "median_better"
        else:
            bayesianVsMedian = "tie"

        counts[f"median_{medianVerdict}"] += 1
        counts[f"bayesian_{bayesianVerdict}"] += 1
        counts[bayesianVsMedian] += 1

        writer.writerow([
            row["EdgeId"], row["StepIndex"], anchorIdA, anchorIdB,
            f"ok:{len(bestList)}/{len(intervals)}:{referenceName}", len(truthSubstring),
            distanceNoRepair, distanceMedian, distanceBayesian,
            medianVerdict, bayesianVerdict, bayesianVsMedian,
            truthSubstring, consensusNoRepair, consensusMedian, consensusBayesian])

print()
print("Against no-repair baseline:")
print("  Median:   helped =", counts["median_helped"], " hurt =", counts["median_hurt"],
    " neutral =", counts["median_neutral"])
print("  Bayesian: helped =", counts["bayesian_helped"], " hurt =", counts["bayesian_hurt"],
    " neutral =", counts["bayesian_neutral"])
print()
print("Bayesian vs Median directly:")
print("  bayesian_better:", counts["bayesian_better"])
print("  median_better:  ", counts["median_better"])
print("  tie:            ", counts["tie"])
print("  skipped:        ", counts["skipped"] + counts["skipped_too_long"])
print()
print("Report written to", arguments.output)
