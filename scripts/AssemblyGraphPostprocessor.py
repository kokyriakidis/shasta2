#!/usr/bin/python3

from shasta2 import *

"""
This demonstrates scripting on the AssemblyGraph.
This must run with the current directory set to the assembly directory
and requires binary data to be available.
"""

# Open the file to get performance log information.
openPerformanceLog("Python-performance.log")

# Get the options.
options = Options()

# Create the Assembler and access what we need.
assembler = Assembler()
assembler.accessMarkers()
assembler.accessAnchors()
assembler.accessJourneys()

# Create a Detangler. Only needed if we want to perform
# detangling operations.
detangler = LikelihoodRatioDetangler(
    options.detangleMinCommonCoverage,
    options.detangleEpsilon,
    options.detangleMaxLogP,
    options.detangleMinLogPDelta,
    options.detangleHighCoverageThreshold)
    
# Read our AssemblyGraph to work on.
assemblyGraph = assembler.getAssemblyGraph("Z", options)

# Here we can operate on this assembly graph.

# When done, we can write the new AssemblyGraph.
assemblyGraph.write("ScriptOutput")
