#!/usr/bin/python3

import shasta2


# Get the arguments.
import argparse
parser = argparse.ArgumentParser(description = "Run a single one detangling iteration.")
parser.add_argument("inputStage", type=str, help="Input assembly stage.")
parser.add_argument("outputStage", type=str, help="Output assembly stage.")
arguments = parser.parse_args()

shasta2.openPerformanceLog("Python-performance.log")

# Get the options from shasta2.conf.
options = shasta2.Options()

assembler = shasta2.Assembler()
assembler.accessAnchors()
assembler.accessJourneys()

# Load the specified assembly stage and assemble sequence.
assemblyGraph = assembler.getAssemblyGraph(arguments.inputStage, options)
assemblyGraph.detangleVertices()
assemblyGraph.strandSymmetricCompress()
assemblyGraph.write(arguments.outputStage)





