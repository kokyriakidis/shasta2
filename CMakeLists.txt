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



# Write out what we will actually build.
message(STATUS "BUILD_EXECUTABLE is " ${BUILD_EXECUTABLE})
message(STATUS "BUILD_PYTHON_MODULE is " ${BUILD_PYTHON_MODULE})



# Option to build with -march=native.
option(BUILD_NATIVE "Build with -march=native." OFF)
message(STATUS "BUILD_NATIVE is " ${BUILD_NATIVE})



# Option to request a debug build.
option(BUILD_DEBUG "Make a debuggable build." OFF)
message(STATUS "BUILD_DEBUG is " ${BUILD_DEBUG})



# The BUILD_ID can be specified to identify the build
# This is normally used only when building a new GitHub release,
# in which case we use the following option when running Cmake:
# -DBUILD_ID="Shasta Release X.Y.Z"
# BUILD_ID should not contain a comma, otherwise the C preprocessor does not handle it correctly.
if(NOT DEFINED BUILD_ID)
    set(BUILD_ID "Shasta development build. This is not a released version.")
endif(NOT DEFINED BUILD_ID)
message(STATUS "BUILD_ID is: " ${BUILD_ID})


# Option to request a build with long markers.
option(SHASTA_LONG_MARKERS "Make a build with SHASTA_LONG_MARKERS set." OFF)
message(STATUS "SHASTA_LONG_MARKERS is " ${SHASTA_LONG_MARKERS})



# Add the subdirectories we need.
if(BUILD_EXECUTABLE)
    add_subdirectory(Executable)
endif(BUILD_EXECUTABLE)

if(BUILD_PYTHON_MODULE)
    add_subdirectory(PythonModule)
endif(BUILD_PYTHON_MODULE)



# Install to the shasta-install directory.
set(CMAKE_INSTALL_PREFIX .)

# Install the scripts.
file(GLOB SCRIPTS scripts/*.py scripts/*.sh)
install(PROGRAMS ${SCRIPTS} DESTINATION shasta-install/bin)





