#!/usr/bin/python3

from shasta2 import *

# Get the arguments.
import argparse
parser = argparse.ArgumentParser()
# parser.add_argument('edgeId0', type=int)
# parser.add_argument('direction', type=int)
parser.add_argument('lowCoverageThreshod', type=int)
parser.add_argument('highCoverageThreshod', type=int)
arguments = parser.parse_args()

openPerformanceLog("Python-performance.log")
options = Options()

# Create the Assembler and access what we need.
assembler = Assembler()
assembler.accessMarkers()
assembler.accessAnchors()
assembler.accessJourneys()


assemblyGraph = assembler.getAssemblyGraph("Z", options)
# assemblyGraph.testLocalSearch(arguments.edgeId0, arguments.direction, arguments.lowCoverageThreshod, arguments.highCoverageThreshod);
assemblyGraph.createSearchGraph(arguments.lowCoverageThreshod, arguments.highCoverageThreshod)
assemblyGraph.assembleAll(0);
assemblyGraph.write("X");
assemblyGraph.writeFasta("X");



