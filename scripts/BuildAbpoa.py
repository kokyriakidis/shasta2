#!/usr/bin/python3

"""
This clones the abpoa repository and builds the abpoa static
and shared libraries. 

Debian/Ubuntu package abpoa only provides the abpoa executable.
abpoa releases also only provide the abpoa executable.

Output goes to ~/.shasta2Build/abpoa
Include files in ~/.shasta2Build/abpoa/abPOA/include
Libraries in ~/.shasta2Build/abpoa/abPOA/lib/

The abpoa repository includes both a CMakeLists.txt and a Makefile.
However we can't use either because:
* CMakeLists.txt uses an unconditional -march=native which 
  results in a non-portable library.
* Makefile does not support building a shared library.  

"""

import os

# Get the path to the root directory.
rootDirectoryWithTilda = "~/.shasta2Build"
rootDirectory = os.path.expanduser(rootDirectoryWithTilda)

# Create the root directory, if it does not already exists.
if os.path.exists(rootDirectory):
    if not os.path.isdir(rootDirectory):
        raise Exception(rootDirectory + " exists but is not a directory")
else:
    os.mkdir(rootDirectory)
    
# cd to the root directory.
os.chdir(rootDirectory)

# Create the abpoa directory, if it does not already exist.
if os.path.exists("abpoa"):
    raise Exception(rootDirectory + "/" + "abpoa already exists.")
else:
    os.mkdir("abpoa")

# cd to the abpoa directory.
os.chdir("abpoa")

# Clone the abpoa repository.
os.system("git clone https://github.com/yangao07/abPOA.git")

# cd to the abPOA directory.
os.chdir("abPOA")

# Create the lib directory.
os.mkdir("lib")

# Change to the src directory.
os.chdir("src")

# Common flags for all compilations.
commonFlags = "-I ../include -O3 -Wall -Wno-unused-function -Wno-misleading-indentation -Wno-stringop-overflow -fno-tree-vectorize"

# Build the static library.
print("Building the abpoa static library.")
os.system("cc -c " + commonFlags + " *.c")
os.system("ar -csr ../lib/libapoa.a *.o")
os.system("rm *.o")

# Build the shared library
print("Building the abpoa shared library.")
os.system("cc -c -fPIC " + commonFlags + " *.c")
os.system("cc -shared -o ../lib/libapoa.so *.o")
os.system("rm *.o")

