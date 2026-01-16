#!/usr/bin/python3

# This installs dependencies required to build shasta2.

import os

# Install the packages we need.
packages = [
    "cmake",
    "g++",
    "libboost-all-dev",
    "libcli11-dev",
    "libpng-dev",
    "python3-dev",
    "python3-pybind11",
    "pybind11-dev",
    ]    
command = "sudo apt-get install " + " ".join(packages)
os.system(command)

# Build abpoa and poasta
import BuildAbpoa
import BuildPoasta


