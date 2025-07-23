#!/usr/bin/python3

import shasta2

options = shasta2.Options()
shasta2.openPerformanceLog("Python-performance.log")
assembler = shasta2.Assembler()
assembler.accessMarkers()
assembler.createMarkerKmers(options.maxMarkerErrorRate, 0)

