#!/usr/bin/python3

"""
This  installs the Rust build environment (in ~.cargo, ~.rustup, then 
clones the kokyriakidis/poasta-c repository and builds 
the poasta static and shared libraries. 

Output goes to ~/.shasta2Build/poasta-c

"""

import os

# Install the rust build environment.
print("Instaling the rust build environment in ~/.cargo, ~/.rustup.")
os.system("curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y")
os.environ['PATH'] = os.environ['PATH'] + ":" + os.path.expanduser("~") + "/.cargo/bin"

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

# Create the poasta-c directory, if it does not already exist.
if os.path.exists("poasta-c"):
    raise Exception(rootDirectory + "/" + "poasta-c already exists.")
else:
    os.mkdir("poasta-c")

# cd to the poasta-c directory.
os.chdir("poasta-c")

# Clone the kokyriakidis/poasta-c repository.
os.system("git clone https://github.com/kokyriakidis/poasta-c.git")

# cd to the poasta-c directory.
os.chdir("poasta-c")

# Build
os.system("cargo build --release")


