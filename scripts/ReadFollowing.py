#!/usr/bin/python3

"""
This uses read following in the complete AnchorGraph
to create the AnchorGraph to be used for assembly.
"""

import shasta2

shasta2.openPerformanceLog("Python-performance.log")

options = shasta2.Options()

assembler = shasta2.Assembler()
assembler.accessMarkers()
assembler.accessMarkerKmers()
assembler.accessAnchors()
assembler.accessJourneys()
assembler.accessCompleteAnchorGraph()
assembler.readFollowing()
