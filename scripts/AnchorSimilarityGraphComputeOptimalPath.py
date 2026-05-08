#!/usr/bin/python3

import shasta2

# Get the AnchorId.
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('componentId', type=int)
arguments = parser.parse_args()

# Access what we need.
assembler = shasta2.Assembler()
assembler.accessMarkers()
assembler.accessMarkerKmers()
assembler.accessAnchors()
assembler.accessAnchorSimilarityGraph()

# Call the shortest path code.
assembler.anchorSimilarityGraphComputeOptimalPath(arguments.componentId)