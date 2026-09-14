#!/usr/bin/python3

# Part of the msa1 hard-region evaluation harness. Takes the outputs of
# FindMsa1HardRegions.py and EvaluateMsa1AgainstTruth.py and checks, for each
# candidate region, whether hifiasm's own assembly (run separately on the same
# reads, see the workflow this script's --hifiasm-fasta expects) already has
# the answer msa1 is trying to vote its way to - by mapping the same core
# reads directly onto hifiasm's contigs (instead of onto HG002 v1.1 truth),
# lifting anchorIdA/anchorIdB over the same way, and comparing hifiasm's
# substring against the "Truth" column EvaluateMsa1AgainstTruth.py already
# established. No dependency on HG002 v1.1 coordinates: hifiasm's sequence is
# compared directly to the truth substring already on hand.

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
    if len(a) * len(b) > cap:
        return None
    previousRow = list(range(len(b) + 1))
    for i, ca in enumerate(a, 1):
        currentRow = [i] + [0] * len(b)
        for j, cb in enumerate(b, 1):
            currentRow[j] = min(
                previousRow[j] + 1,
                currentRow[j - 1] + 1,
                previousRow[j - 1] + (ca != cb))
        previousRow = currentRow
    return previousRow[-1]

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
    "Check whether hifiasm's own assembly (mapped the same way as the HG002 "
    "truth) already agrees with the answer msa1 votes for at each candidate "
    "region.")
parser.add_argument("stage", type=str)
parser.add_argument("--candidates", type=str, default="Msa1CandidateRegions.csv")
parser.add_argument("--truth-report", type=str, default="Msa1TruthReport.csv")
parser.add_argument("--hifiasm-fasta", type=str, required=True,
    help="FASTA of hifiasm's own contig(s) - e.g. both haplotype p_ctg files "
         "concatenated - built by running hifiasm separately on the same reads.")
parser.add_argument("--minimap2", type=str,
    default=os.path.expanduser("~/Downloads/minimap2-2.30_x64-linux/minimap2"))
parser.add_argument("--output", type=str, default="Msa1VsHifiasm.csv")
parser.add_argument("--work-dir", type=str, default=None)
arguments = parser.parse_args()

options = shasta2.Options()
assembler = shasta2.Assembler()
assembler.accessAnchors()
assembler.accessJourneys()

candidates = {}
with open(arguments.candidates, newline="") as csvFile:
    for row in csv.DictReader(csvFile):
        candidates[(row["EdgeId"], row["StepIndex"])] = row

rows = []
with open(arguments.truth_report, newline="") as csvFile:
    for row in csv.DictReader(csvFile):
        if row["TruthStatus"].startswith("ok") and row["DistanceAdaptive"] != "":
            rows.append(row)
print(len(rows), "regions with established truth loaded from", arguments.truth_report)

for row in rows:
    anchorIdA = int(row["AnchorIdA"])
    anchorIdB = int(row["AnchorIdB"])
    orientedReadIds = candidates[(row["EdgeId"], row["StepIndex"])]["OrientedReadIds"].split(";")
    row["CoreReads"] = [
        orientedReadId for orientedReadId in orientedReadIds
        if assembler.anchorContainsOrientedRead(anchorIdA, orientedReadId)
        and assembler.anchorContainsOrientedRead(anchorIdB, orientedReadId)]

allReads = sorted({orientedReadId for row in rows for orientedReadId in row["CoreReads"]})
print(len(allReads), "distinct oriented reads to map onto hifiasm's assembly.")

workDir = arguments.work_dir or tempfile.mkdtemp(prefix="msa1VsHifiasm_")
os.makedirs(workDir, exist_ok=True)

readSequences = {}
readsFastaName = os.path.join(workDir, "reads.fasta")
with open(readsFastaName, "w") as fasta:
    for orientedReadId in allReads:
        sequence = assembler.getOrientedReadSequenceString(orientedReadId)
        readSequences[orientedReadId] = sequence
        fasta.write(f">{orientedReadId}\n{sequence}\n")

samFileName = os.path.join(workDir, "reads.sam")
with open(samFileName, "w") as samFile:
    subprocess.run(
        [arguments.minimap2, "-a", "--eqx", "-x", "map-ont",
            arguments.hifiasm_fasta, readsFastaName],
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
        alignments[queryName] = (fields[2], int(fields[3]) - 1, fields[5], bool(flag & 16))
print(len(alignments), "of", len(allReads), "reads have a primary mapping onto hifiasm's assembly.")

with open(arguments.hifiasm_fasta) as f:
    contigSequences = {}
    name = None
    chunks = []
    for line in f:
        line = line.rstrip("\n")
        if line.startswith(">"):
            if name is not None:
                contigSequences[name] = "".join(chunks)
            name = line[1:].split()[0]
            chunks = []
        else:
            chunks.append(line)
    if name is not None:
        contigSequences[name] = "".join(chunks)



fixedCount = hifiasmAgreesAtFixed = 0
brokeCount = hifiasmAgreesAtBroke = 0
hifiasmExactCount = 0

with open(arguments.output, "w", newline="") as outputFile:
    writer = csv.writer(outputFile)
    writer.writerow([
        "EdgeId", "StepIndex", "DistanceAdaptive", "DistanceMsa1",
        "LengthAdaptive", "LengthMsa1",
        "HifiasmStatus", "DistanceHifiasm", "HifiasmLength", "Category"])

    for row in rows:
        anchorIdA = int(row["AnchorIdA"])
        anchorIdB = int(row["AnchorIdB"])
        distanceAdaptive = int(row["DistanceAdaptive"])
        distanceMsa1 = int(row["DistanceMsa1"])
        lengthAdaptive = int(row["LengthAdaptive"])
        lengthMsa1 = int(row["LengthMsa1"])
        truth = row["Truth"]

        category = (
            "fixed" if distanceAdaptive > 0 and distanceMsa1 == 0 else
            "broke" if distanceAdaptive == 0 and distanceMsa1 > 0 else
            "other")
        if category == "fixed":
            fixedCount += 1
        elif category == "broke":
            brokeCount += 1

        intervals = []
        for orientedReadId in row["CoreReads"]:
            alignment = alignments.get(orientedReadId)
            if alignment is None:
                continue
            contigName, contigStart, cigar, isReverse = alignment

            positionA = assembler.getAnchorPositionInOrientedRead(anchorIdA, orientedReadId)
            positionB = assembler.getAnchorPositionInOrientedRead(anchorIdB, orientedReadId)
            sequenceLength = len(readSequences[orientedReadId])

            queryPositionA = (sequenceLength - 1 - positionA) if isReverse else positionA
            queryPositionB = (sequenceLength - 1 - positionB) if isReverse else positionB

            contigPositionA = queryToReferencePosition(cigar, contigStart, queryPositionA)
            contigPositionB = queryToReferencePosition(cigar, contigStart, queryPositionB)
            if contigPositionA is None or contigPositionB is None:
                continue

            if isReverse:
                begin = min(contigPositionA, contigPositionB) + 1
                end = max(contigPositionA, contigPositionB) + 1
            else:
                begin = min(contigPositionA, contigPositionB)
                end = max(contigPositionA, contigPositionB)
            if begin >= end:
                continue
            intervals.append((contigName, begin, end, isReverse))

        if not intervals:
            writer.writerow([row["EdgeId"], row["StepIndex"], distanceAdaptive, distanceMsa1,
                lengthAdaptive, lengthMsa1, "no_reads_mapped", "", "", category])
            continue

        groups = defaultdict(list)
        for contigName, begin, end, isReverse in intervals:
            groups[(contigName, isReverse)].append((begin, end))
        (contigName, isReverse), bestList = max(groups.items(), key=lambda kv: len(kv[1]))
        if len(bestList) < max(1, len(intervals) // 2 + 1):
            writer.writerow([row["EdgeId"], row["StepIndex"], distanceAdaptive, distanceMsa1,
                lengthAdaptive, lengthMsa1,
                f"reads_disagree:{len(bestList)}/{len(intervals)}", "", "", category])
            continue

        begins = sorted(b for b, e in bestList)
        ends = sorted(e for b, e in bestList)
        begin = begins[len(begins) // 2]
        end = ends[len(ends) // 2]

        hifiasmSlice = contigSequences[contigName][begin:end]
        if isReverse:
            hifiasmSlice = reverseComplement(hifiasmSlice)

        row["HifiasmLength"] = len(hifiasmSlice)
        distanceHifiasm = editDistance(hifiasmSlice, truth)
        if distanceHifiasm == 0:
            hifiasmExactCount += 1
        # "Correct" means exactly matching truth - for fixed regions that's
        # msa1's own answer, for broke regions that's adaptive's.
        if category == "fixed" and distanceHifiasm == 0:
            hifiasmAgreesAtFixed += 1
        if category == "broke" and distanceHifiasm == 0:
            hifiasmAgreesAtBroke += 1

        writer.writerow([
            row["EdgeId"], row["StepIndex"], distanceAdaptive, distanceMsa1,
            lengthAdaptive, lengthMsa1,
            f"ok:{len(bestList)}/{len(intervals)}:{contigName}",
            distanceHifiasm, row["HifiasmLength"], category])

print()
print("fixed regions (msa1 exact, adaptive wrong):", fixedCount,
      " hifiasm also exact:", hifiasmAgreesAtFixed)
print("broke regions (adaptive exact, msa1 wrong):", brokeCount,
      " hifiasm also exact:", hifiasmAgreesAtBroke)
print("hifiasm exact matches to truth (all categories):", hifiasmExactCount, "/", len(rows))
print("Report written to", arguments.output)
