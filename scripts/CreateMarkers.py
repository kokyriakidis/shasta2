#!/usr/bin/python3

import shasta2

shasta2.openPerformanceLog("Python-performance.log")
options = shasta2.Options()
assembler = shasta2.Assembler()
assembler.createKmerChecker(options.k, options.markerDensity)
assembler.createMarkers()

