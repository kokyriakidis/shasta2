#!/usr/bin/python3

import shasta2

# Get the AnchorId.
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('anchorId', type=str)
arguments = parser.parse_args()
anchorIdString = arguments.anchorId
anchorId = shasta2.anchorIdFromString(anchorIdString);

# Access what we need.
assembler = shasta2.Assembler()
assembler.accessMarkers()
assembler.accessMarkerKmers()
assembler.accessAnchors()
assembler.accessAnchorSimilarityGraph()

# Call the shortest path code.
assembler.anchorSimilarityGraphCreateShortestPathTree(anchorId)