# For building, the following are required:
# - The Boost C++ libraries (Ubuntu package libboost-dev).
# - Cli11 (Ubuntu package libcli11-dev).
# - The Abpoa library. Ubuntu package abpoa only provides the executable.
#   See ~/abpoa/README for build directions.
#   * The include files must be in ~/abpoa/abPOA/include. 
#   * The static and dynamic libraries (libabpoa.a and libabpoa.so) must be in /usr/local/lib.
# For running, the following are required:
# - Poasta executable must be in the PATH. See ~/poasta/README for build directions.
# - If using the http server:
#   * Graphviz (Ubuntu package graphviz).
#   * ogdf-fmmm (see ~/ogdf-fmmm/ogdf-fmm/README for build directions). 
#     Must be on path under name customLayout.


cmake_minimum_required(VERSION 3.16)
project(shasta2)


# Figure out if we are on Intel (as opposed to ARM).
set(X86_64 OFF)
if(${CMAKE_SYSTEM_PROCESSOR} MATCHES "x86_64")
    set(X86_64 ON)
endif()



# Summarize the environment.
message(STATUS "Processor architecture is " ${CMAKE_HOST_SYSTEM_PROCESSOR})
message(STATUS "CMAKE_SYSTEM_NAME is " ${CMAKE_SYSTEM_NAME})



# Decide what we want to build.
option(BUILD_EXECUTABLE "Build the executable." ON)
option(BUILD_PYTHON_MODULE "Build the Python module." ON)
option(BUILD_STATIC_LIBRARY "Build the static library." OFF)



# Write out what we will actually build.
message(STATUS "BUILD_EXECUTABLE is " ${BUILD_EXECUTABLE})
message(STATUS "BUILD_PYTHON_MODULE is " ${BUILD_PYTHON_MODULE})
message(STATUS "BUILD_STATIC_LIBRARY is " ${BUILD_STATIC_LIBRARY})



# Option to request a debug build.
option(BUILD_DEBUG "Make a debuggable build." OFF)
message(STATUS "BUILD_DEBUG is " ${BUILD_DEBUG})



# Add the subdirectories we need.
if(BUILD_EXECUTABLE)
    add_subdirectory(Executable)
endif(BUILD_EXECUTABLE)

if(BUILD_PYTHON_MODULE)
    add_subdirectory(PythonModule)
endif(BUILD_PYTHON_MODULE)

if(BUILD_STATIC_LIBRARY)
    add_subdirectory(StaticLibrary)
endif(BUILD_STATIC_LIBRARY)



# Install to the shasta-install directory.
set(CMAKE_INSTALL_PREFIX .)

# Install the scripts.
file(GLOB SCRIPTS scripts/*.py scripts/*.sh)
install(PROGRAMS ${SCRIPTS} DESTINATION shasta2-install/bin)

# Install the configuration files.
install(DIRECTORY conf DESTINATION shasta2-install)





