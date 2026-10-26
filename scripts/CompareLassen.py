#!/usr/bin/python3

# Compares lassen runs (paoloshasta/lassen) of several reassemblies of the same
# Final assembly graph (AssembleFinalWithMethod.py), as in commit d09909e6:
# total errors (mismatch + insert + delete) over the 25 longest assembled
# segments that are at least 100 kb and map exactly once (one haplotype, not
# fragmented: lassen's Summary2.csv, Number of mappings 1), the same segments
# in every run, and over every mapping of every run (Summary1.csv).
#
# usage: CompareLassen.py name=lassenDirectory ... [--count 25] [--min-length 100000]
# The first run chooses the segments; the others are compared against it.

import argparse
import csv
import os



def errors(row):
    return int(row["Mismatch count"]) + int(row["Insert count"]) + int(row["Delete count"])



def readSummary2(directory):
    rows = {}
    with open(os.path.join(directory, "Summary2.csv")) as f:
        for row in csv.DictReader(f):
            rows[row["Assembled name"]] = row
    return rows



def readSummary1(directory):
    with open(os.path.join(directory, "Summary1.csv")) as f:
        return list(csv.DictReader(f))



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("runs", nargs="+")
    parser.add_argument("--count", type=int, default=25)
    parser.add_argument("--min-length", type=int, default=100000)
    arguments = parser.parse_args()
    runs = [r.split("=", 1) for r in arguments.runs]

    segments = {name: readSummary2(directory) for name, directory in runs}
    first = runs[0][0]
    chosen = [s for s, row in segments[first].items()
              if row["Number of mappings"] == "1" and int(row["Assembled length"]) >= arguments.min_length]
    chosen.sort(key=lambda s: -int(segments[first][s]["Assembled length"]))
    chosen = chosen[:arguments.count]
    for name, _ in runs:
        chosen = [s for s in chosen if segments[name].get(s, {}).get("Number of mappings") == "1"]
    print("%d segments (longest %d, shortest %d bases)" % (len(chosen),
        int(segments[first][chosen[0]]["Assembled length"]), int(segments[first][chosen[-1]]["Assembled length"])))

    print("\n%-10s %10s %9s %8s %8s %12s" % ("run", "errors", "mismatch", "insert", "delete", "change"))
    base = None
    for name, _ in runs:
        rows = [segments[name][s] for s in chosen]
        total = sum(errors(r) for r in rows)
        base = total if base is None else base
        print("%-10s %10d %9d %8d %8d %+11.1f%%" % (name, total,
            sum(int(r["Mismatch count"]) for r in rows), sum(int(r["Insert count"]) for r in rows),
            sum(int(r["Delete count"]) for r in rows), 100.0 * (total - base) / max(1, base)))

    for name, _ in runs[1:]:
        better = sum(1 for s in chosen if errors(segments[name][s]) < errors(segments[first][s]))
        worse = sum(1 for s in chosen if errors(segments[name][s]) > errors(segments[first][s]))
        print("%s against %s: better in %d segments, worse in %d, equal in %d" % (name, first, better, worse,
            len(chosen) - better - worse))

    print("\nEvery mapping")
    print("%-10s %9s %14s %10s %9s %8s %8s %12s" % ("run", "mappings", "aligned bases", "errors", "mismatch",
        "insert", "delete", "change"))
    base = None
    for name, directory in runs:
        rows = readSummary1(directory)
        total = sum(errors(r) for r in rows)
        base = total if base is None else base
        print("%-10s %9d %14d %10d %9d %8d %8d %+11.1f%%" % (name, len(rows),
            sum(int(r["Assembled aligned length"]) for r in rows), total,
            sum(int(r["Mismatch count"]) for r in rows), sum(int(r["Insert count"]) for r in rows),
            sum(int(r["Delete count"]) for r in rows), 100.0 * (total - base) / max(1, base)))



main()
