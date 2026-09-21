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
# compared (by edit distance) to both consensusNoRepair and consensusWithRepair
# to see whether the repair moved the consensus closer to or further from truth.

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



# An O(len(a) * len(b)) edit distance, computed in C++ (shasta2.editDistance)
# rather than here: candidate regions are local homopolymer-repair windows and
# individually short, but there can be thousands of them, and CPython's
# per-cell interpreter overhead made this the bottleneck of the whole
# pipeline. If a pair is too large to make the comparison cheap, it is
# skipped rather than slowed down, same as before.
def editDistance(a, b, cap=4_000_000):
    distance = shasta2.editDistance(a, b, cap)
    return None if distance < 0 else distance



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
parser.add_argument("--window-pad", type=int, default=5000,
    help="Bases of read context kept on each side of anchorIdA/anchorIdB when mapping "
         "a core read to truth (see the comment above the mapping step).")
parser.add_argument("--min-mapq", type=int, default=10,
    help="Minimum minimap2 MAPQ for a core read's mapping to be trusted (see the "
         "comment above the mapping step).")
parser.add_argument("--threads", type=int, default=os.cpu_count() or 4,
    help="Threads to give minimap2.")
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
readSequences = {
    orientedReadId: assembler.getOrientedReadSequenceString(orientedReadId)
    for orientedReadId in allReads}
print(len(allReads), "distinct oriented reads among the candidate regions.")



# Map, per region, only a padded window of each core read around that
# region's own anchorIdA/anchorIdB - not the whole read - and map every
# region's windows together in one minimap2 call. A core read can be tens of
# kb long while every region here is a local homopolymer-repair window a few
# hundred bases wide; mapping full reads made minimap2 the bottleneck of the
# whole pipeline (tens of minutes at real-genome candidate counts, dwarfing
# every other step), purely because its cost tracks total input bases and
# full reads made that far bigger than the problem needed. A read used by
# several regions gets one window per region - windows are keyed by
# (row index, orientedReadId), not just orientedReadId, since the same read
# can need a different window in each region it contributes to.
#
# A smaller window is also a less unique one, and an unqualified mapping is
# worse than no mapping at all for a truth harness: measured against full-read
# mapping as ground truth on a real 3000-window sample, a --window-pad of 1000
# put only 51.7% of windows at MAPQ>=10, but among those, agreement with the
# full-read locus was 100% (0/1400) - MAPQ<10 is where essentially all of the
# disagreement lives (49.4% wrong at MAPQ=0, 14.5% wrong at MAPQ 1-9). Below,
# --min-mapq drops anything under that bar from the truth vote, same as an
# unmapped read - silently trusting a low-MAPQ window would occasionally
# extract truth from the wrong genomic locus entirely (one sampled case
# mapped to a different chromosome), corrupting the comparison without any
# sign anything was wrong. The default --window-pad of 5000 exists to reduce
# how much this filter throws away: at that width MAPQ>=10 covers 59.9% of
# windows instead of 51.7%, while the mapping call itself is still a small
# fraction of full-read mapping's cost (total input bases stay two orders of
# magnitude below mapping every core read in full).
workDir = arguments.work_dir or tempfile.mkdtemp(prefix="msa1Eval_")
os.makedirs(workDir, exist_ok=True)

# windowInfo[(rowIndex, orientedReadId)] = (localPositionA, localPositionB, windowLength):
# that read's anchorIdA/anchorIdB positions in this row, re-based to the
# window's own coordinates, for use once the window's alignment is known.
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

# Keep only the primary alignment of each window (skip unmapped/secondary/supplementary).
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
        mapq = int(fields[4])
        cigar = fields[5]
        alignments[queryName] = (referenceName, referenceStart, cigar, bool(flag & 16), mapq)
trustedAlignments = sum(1 for a in alignments.values() if a[4] >= arguments.min_mapq)
print(len(alignments), "of", len(windowInfo), "windows have a primary mapping to truth,",
    trustedAlignments, f"at MAPQ>={arguments.min_mapq} or better.")



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
        "TruthStatus", "TruthLength", "LengthNoRepair", "LengthWithRepair",
        "DistanceNoRepair", "DistanceWithRepair", "Verdict",
        "Truth", "ConsensusNoRepair", "ConsensusWithRepair"])

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
                # Not unmapped, but not trustworthy either - see the comment
                # above the mapping step. Treated the same as unmapped: this
                # read just doesn't get a vote for this region's truth locus.
                continue

            localPositionA, localPositionB, windowLength = windowInfo[(rowIndex, orientedReadId)]

            # The SAM CIGAR walks the query in the orientation minimap2 aligned it in,
            # which is the reverse complement of our own oriented read sequence when
            # the alignment is on the reverse strand.
            queryPositionA = (windowLength - 1 - localPositionA) if isReverse else localPositionA
            queryPositionB = (windowLength - 1 - localPositionB) if isReverse else localPositionB

            referencePositionA = queryToReferencePosition(cigar, referenceStart, queryPositionA)
            referencePositionB = queryToReferencePosition(cigar, referenceStart, queryPositionB)
            if referencePositionA is None or referencePositionB is None:
                continue

            # LocalAssembly7 assembles a core read's window as
            # [positionA, positionB) - inclusive of anchorIdA's position,
            # EXCLUSIVE of anchorIdB's (LocalAssembly7.cpp:377-378,408: the
            # loop is `for(position=positionBegin; position!=positionEnd; ...)`
            # with positionBegin=positionA, positionEnd=positionB). Reproduce
            # that exactly, strand-aware: forward, referencePositionA <
            # referencePositionB and we want [refA, refB); reverse, the read's
            # own increasing direction runs backward in reference coordinates,
            # so referencePositionB < referencePositionA and the equivalent
            # forward-strand slice (before its later reverse-complementing) is
            # [refB + 1, refA + 1) - inclusive of refA, exclusive of refB.
            if isReverse:
                begin = min(referencePositionA, referencePositionB) + 1
                end = max(referencePositionA, referencePositionB) + 1
            else:
                begin = min(referencePositionA, referencePositionB)
                end = max(referencePositionA, referencePositionB)
            if begin >= end:
                continue
            intervals.append((referenceName, begin, end, isReverse))

        if not intervals:
            writer.writerow([
                row["EdgeId"], row["StepIndex"], anchorIdA, anchorIdB,
                "no_reads_mapped", "", "", "", "", "", "skipped",
                "", row["ConsensusNoRepair"], row["ConsensusWithRepair"]])
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
                "", row["ConsensusNoRepair"], row["ConsensusWithRepair"]])
            skippedCount += 1
            continue

        begins = sorted(b for b, e in bestList)
        ends = sorted(e for b, e in bestList)
        begin = begins[len(begins) // 2]
        end = ends[len(ends) // 2]

        truthSubstring = fetchTruthSubstring(referenceName, begin, end)
        if isReverse:
            truthSubstring = reverseComplement(truthSubstring)

        consensusNoRepair = row["ConsensusNoRepair"]
        consensusWithRepair = row["ConsensusWithRepair"]
        distanceNoRepair = editDistance(consensusNoRepair, truthSubstring)
        distanceWithRepair = editDistance(consensusWithRepair, truthSubstring)

        if distanceNoRepair is None or distanceWithRepair is None:
            verdict = "skipped_too_long"
            skippedCount += 1
        elif distanceWithRepair < distanceNoRepair:
            verdict = "helped"
            helpedCount += 1
        elif distanceWithRepair > distanceNoRepair:
            verdict = "hurt"
            hurtCount += 1
        else:
            verdict = "neutral"
            neutralCount += 1

        writer.writerow([
            row["EdgeId"], row["StepIndex"], anchorIdA, anchorIdB,
            f"ok:{len(bestList)}/{len(intervals)}:{referenceName}", len(truthSubstring),
            len(consensusNoRepair), len(consensusWithRepair),
            distanceNoRepair, distanceWithRepair, verdict,
            truthSubstring, consensusNoRepair, consensusWithRepair])

print()
print("Verdict summary:")
print("  helped: ", helpedCount)
print("  hurt:   ", hurtCount)
print("  neutral:", neutralCount)
print("  skipped:", skippedCount)
print("Report written to", arguments.output)
