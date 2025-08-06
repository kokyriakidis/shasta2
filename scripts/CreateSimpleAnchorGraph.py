#!/usr/bin/python3

import shasta2

shasta2.openPerformanceLog("Python-performance.log")

assembler = shasta2.Assembler()
assembler.accessMarkers()
assembler.accessAnchors()
assembler.accessJourneys()
assembler.createSimpleAnchorGraph()
