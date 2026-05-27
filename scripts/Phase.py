#!/usr/bin/python3

import shasta2

# Get the argument.
import argparse
parser = argparse.ArgumentParser(description = "Do a round oh phasing on a given assembly stage.")
parser.add_argument("stage", type=str, help="Assembly stage.")
arguments = parser.parse_args()

shasta2.openPerformanceLog("Python-performance.log")

# Get the options from shasta2.conf.
options = shasta2.Options()

# Create the Assembler and access what we need.
assembler = shasta2.Assembler()
assembler.accessAnchors()
assembler.accessJourneys()

# Load the specified assembly stage and assemble sequence.
assemblyGraph = assembler.getAssemblyGraph(arguments.stage, options)
assemblyGraph.phaseSuperbubbleChains()

# Write it out.
assembledName = arguments.stage + "-Phased"
assemblyGraph.write(assembledName)




