#!/usr/bin/python3

# Part of the msa1 hard-region evaluation harness. See FindMsa1HardRegions.py,
# which produces the --candidates csv this script reads.
#
# For every candidate region (msa1's repair changed the consensus), establish
# ground truth by mapping every "core" contributing read - every oriented read
# that is on both anchorIdA and anchorIdB, i.e. every read AnchorPair itself
# would have used, as opposed to the extra context reads
# AssemblyGraph.getAssemblyGraphSteps() borrows from the previous/next step -
# to the HG002 v1.1 diploid truth assembly with minimap2, then lifting each
# read's anchorIdA/anchorIdB marker positions over to truth-reference
# coordinates by walking that read's CIGAR. Reads that agree on the truth
# interval give the ground truth substring for the region, which is then
# compared (by edit distance) to both consensusAdaptive and consensusMsa1 to
# see whether the repair moved the consensus closer to or further from truth.

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



# A simple O(len(a) * len(b)) edit distance. Candidate regions are local
# homopolymer-repair windows and are expected to be short; if a pair is too
# large to make this cheap, the comparison is skipped rather than slowed down.
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



# CIGAR ops that consume the query (the SAM SEQ field) and/or the reference.
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



# alignmentStart: 0-based reference position of the first base described by the CIGAR.
# queryPosition: 0-based position in the SAM SEQ field (the query as aligned, i.e.
# reverse complemented relative to our own oriented read sequence when the
# alignment is on the reverse strand).
# Returns None if queryPosition falls inside an insertion/soft-clip, where there
# is no single corresponding reference position.
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
    "Compare msa1-repaired consensus against HG002 v1.1 truth for the candidate "
    "regions found by FindMsa1HardRegions.py.")
parser.add_argument("stage", type=str,
    help="Assembly stage the candidates were found in (must be reachable with the "
         "same Data/ directory used by FindMsa1HardRegions.py).")
parser.add_argument("--candidates", type=str, default="Msa1CandidateRegions.csv")
parser.add_argument("--truth-fasta", type=str,
    default=os.path.expanduser("~/Downloads/hg002v1.1.fasta"))
parser.add_argument("--minimap2", type=str,
    default=os.path.expanduser("~/Downloads/minimap2-2.30_x64-linux/minimap2"))
parser.add_argument("--samtools", type=str, default="samtools")
parser.add_argument("--truth-index", type=str, default=None,
    help="Prebuilt minimap2 map-ont index (.mmi) for --truth-fasta. Built once next to "
         "the fasta and reused on later runs by default, so repeated runs don't "
         "re-index the (large) truth genome from scratch every time.")
parser.add_argument("--output", type=str, default="Msa1TruthReport.csv")
parser.add_argument("--work-dir", type=str, default=None,
    help="Directory for the temporary reads fasta/sam files (default: a fresh temp dir).")
arguments = parser.parse_args()

options = shasta2.Options()
assembler = shasta2.Assembler()
assembler.accessAnchors()
assembler.accessJourneys()

# Build (once) or reuse a minimap2 index for the truth fasta, so repeated runs
# of this script don't re-index the whole genome every time.
truthIndex = arguments.truth_index or (arguments.truth_fasta + ".map-ont.mmi")
if not os.path.exists(truthIndex):
    print("Building minimap2 index", truthIndex, "(one-time cost for this truth fasta)...")
    subprocess.run(
        [arguments.minimap2, "-x", "map-ont", "-d", truthIndex, arguments.truth_fasta],
        check=True)



# Load the candidates.
rows = []
with open(arguments.candidates, newline="") as csvFile:
    for row in csv.DictReader(csvFile):
        row["OrientedReadIds"] = row["OrientedReadIds"].split(";") if row["OrientedReadIds"] else []
        rows.append(row)
print(len(rows), "candidate regions loaded from", arguments.candidates)



# For each region, keep only the "core" reads - the ones on both anchorIdA and
# anchorIdB, i.e. the ones AnchorPair itself defines - and not the extra
# context reads borrowed from the previous/next step, which are not
# guaranteed to be on this step's anchors at all.
for row in rows:
    anchorIdA = int(row["AnchorIdA"])
    anchorIdB = int(row["AnchorIdB"])
    row["CoreReads"] = [
        orientedReadId for orientedReadId in row["OrientedReadIds"]
        if assembler.anchorContainsOrientedRead(anchorIdA, orientedReadId)
        and assembler.anchorContainsOrientedRead(anchorIdB, orientedReadId)]

allReads = sorted({orientedReadId for row in rows for orientedReadId in row["CoreReads"]})
print(len(allReads), "distinct oriented reads to map to truth.")



# Write them out and map them all to truth in one minimap2 call.
workDir = arguments.work_dir or tempfile.mkdtemp(prefix="msa1Eval_")
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
            truthIndex, readsFastaName],
        stdout=samFile, stderr=subprocess.DEVNULL, check=True)

# Keep only the primary alignment of each read (skip unmapped/secondary/supplementary).
alignments = {}
with open(samFileName) as samFile:
    for line in samFile:
        if line.startswith("@"):
            continue
        fields = line.rstrip("\n").split("\t")
        queryName = fields[0]
        flag = int(fields[1])
        if flag & 0x904:   # unmapped(0x4) | secondary(0x100) | supplementary(0x800)
            continue
        referenceName = fields[2]
        referenceStart = int(fields[3]) - 1   # SAM POS is 1-based.
        cigar = fields[5]
        alignments[queryName] = (referenceName, referenceStart, cigar, bool(flag & 16))
print(len(alignments), "of", len(allReads), "reads have a primary mapping to truth.")



def fetchTruthSubstring(referenceName, begin, end):
    # begin, end: 0-based, end-exclusive.
    region = f"{referenceName}:{begin + 1}-{end}"
    result = subprocess.run(
        [arguments.samtools, "faidx", arguments.truth_fasta, region],
        capture_output=True, text=True, check=True)
    return "".join(result.stdout.splitlines()[1:])   # Skip the ">" header line.



# For each region, lift each core read's anchorIdA/anchorIdB positions over to
# truth coordinates, cross-check agreement across reads, and extract the
# truth substring.
helpedCount = 0
hurtCount = 0
neutralCount = 0
skippedCount = 0

with open(arguments.output, "w", newline="") as outputFile:
    writer = csv.writer(outputFile)
    writer.writerow([
        "EdgeId", "StepIndex", "AnchorIdA", "AnchorIdB",
        "TruthStatus", "TruthLength", "LengthAdaptive", "LengthMsa1",
        "DistanceAdaptive", "DistanceMsa1", "Verdict",
        "Truth", "ConsensusAdaptive", "ConsensusMsa1"])

    for row in rows:
        anchorIdA = int(row["AnchorIdA"])
        anchorIdB = int(row["AnchorIdB"])

        intervals = []
        for orientedReadId in row["CoreReads"]:
            alignment = alignments.get(orientedReadId)
            if alignment is None:
                continue
            referenceName, referenceStart, cigar, isReverse = alignment

            positionA = assembler.getAnchorPositionInOrientedRead(anchorIdA, orientedReadId)
            positionB = assembler.getAnchorPositionInOrientedRead(anchorIdB, orientedReadId)
            sequenceLength = len(readSequences[orientedReadId])

            # The SAM CIGAR walks the query in the orientation minimap2 aligned it in,
            # which is the reverse complement of our own oriented read sequence when
            # the alignment is on the reverse strand.
            queryPositionA = (sequenceLength - 1 - positionA) if isReverse else positionA
            queryPositionB = (sequenceLength - 1 - positionB) if isReverse else positionB

            referencePositionA = queryToReferencePosition(cigar, referenceStart, queryPositionA)
            referencePositionB = queryToReferencePosition(cigar, referenceStart, queryPositionB)
            if referencePositionA is None or referencePositionB is None:
                continue

            begin = min(referencePositionA, referencePositionB)
            end = max(referencePositionA, referencePositionB) + 1
            intervals.append((referenceName, begin, end, isReverse))

        if not intervals:
            writer.writerow([
                row["EdgeId"], row["StepIndex"], anchorIdA, anchorIdB,
                "no_reads_mapped", "", "", "", "", "", "skipped",
                "", row["ConsensusAdaptive"], row["ConsensusMsa1"]])
            skippedCount += 1
            continue

        # Group by (referenceName, strand) and require a majority to agree.
        groups = defaultdict(list)
        for referenceName, begin, end, isReverse in intervals:
            groups[(referenceName, isReverse)].append((begin, end))
        (referenceName, isReverse), bestList = max(groups.items(), key=lambda kv: len(kv[1]))

        if len(bestList) < max(1, len(intervals) // 2 + 1):
            writer.writerow([
                row["EdgeId"], row["StepIndex"], anchorIdA, anchorIdB,
                f"reads_disagree:{len(bestList)}/{len(intervals)}",
                "", "", "", "", "", "skipped",
                "", row["ConsensusAdaptive"], row["ConsensusMsa1"]])
            skippedCount += 1
            continue

        begins = sorted(b for b, e in bestList)
        ends = sorted(e for b, e in bestList)
        begin = begins[len(begins) // 2]
        end = ends[len(ends) // 2]

        truthSubstring = fetchTruthSubstring(referenceName, begin, end)
        if isReverse:
            truthSubstring = reverseComplement(truthSubstring)

        consensusAdaptive = row["ConsensusAdaptive"]
        consensusMsa1 = row["ConsensusMsa1"]
        distanceAdaptive = editDistance(consensusAdaptive, truthSubstring)
        distanceMsa1 = editDistance(consensusMsa1, truthSubstring)

        if distanceAdaptive is None or distanceMsa1 is None:
            verdict = "skipped_too_long"
            skippedCount += 1
        elif distanceMsa1 < distanceAdaptive:
            verdict = "helped"
            helpedCount += 1
        elif distanceMsa1 > distanceAdaptive:
            verdict = "hurt"
            hurtCount += 1
        else:
            verdict = "neutral"
            neutralCount += 1

        writer.writerow([
            row["EdgeId"], row["StepIndex"], anchorIdA, anchorIdB,
            f"ok:{len(bestList)}/{len(intervals)}:{referenceName}", len(truthSubstring),
            len(consensusAdaptive), len(consensusMsa1),
            distanceAdaptive, distanceMsa1, verdict,
            truthSubstring, consensusAdaptive, consensusMsa1])

print()
print("Verdict summary:")
print("  helped: ", helpedCount)
print("  hurt:   ", hurtCount)
print("  neutral:", neutralCount)
print("  skipped:", skippedCount)
print("Report written to", arguments.output)
