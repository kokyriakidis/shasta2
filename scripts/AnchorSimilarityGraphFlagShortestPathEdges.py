#!/usr/bin/python3

import shasta2

# Access what we need.
assembler = shasta2.Assembler()
assembler.accessMarkers()
assembler.accessAnchors()
assembler.accessAnchorSimilarityGraph()

# Call the shortest path code.
assembler.anchorSimilarityGraphFlagShortestPathEdges()