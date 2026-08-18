#!/usr/bin/python3

import shasta2

# Get the argument.
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("stage", type=str, help="Assembly stage.")
arguments = parser.parse_args()

options = shasta2.Options()

assembler = shasta2.Assembler()
assembler.accessAnchors()
assembler.accessJourneys()

assemblyGraph = assembler.getAssemblyGraph(arguments.stage, options)
print(assemblyGraph.countIsolatedVertices())



