cmake_minimum_required(VERSION 3.16)
project(shasta2Executable)

# C++ dialect.
add_definitions(-std=c++20)

# Compilation warnings.
add_definitions(-Wall -Wconversion -Wno-unused-result -Wno-trigraphs -Wno-psabi -Wunused-parameter)

# Optimization and debug options.
if(BUILD_DEBUG)
    add_definitions(-ggdb3)
    add_definitions(-O0)
else(BUILD_DEBUG)
    add_definitions(-g0)
    add_definitions(-O3)
    # NDEBUG is required to turn off SeqAn debug code.
    add_definitions(-DNDEBUG)
endif(BUILD_DEBUG)

# Definitions needed to eliminate dependency on the boost system library.
add_definitions(-DBOOST_SYSTEM_NO_DEPRECATED)
add_definitions(-DBOOST_ERROR_CODE_HEADER_ONLY)

# This is needed to avoid some Boost warnings about deprecated messages.
# We canot fix this in Shasta as it is a Boost problem.
# It can be removed if Boost is fixed.
add_definitions(-DBOOST_ALLOW_DEPRECATED_HEADERS)

# Source files
file(GLOB SOURCES ../src/*.cpp)

# Include directory.
include_directories(../src)

# Define our executable.
add_executable(shasta2Executable ${SOURCES})
set_target_properties(shasta2Executable PROPERTIES OUTPUT_NAME "shasta2")

# Request a static executable.
set_target_properties(shasta2Executable PROPERTIES LINK_FLAGS "-static" )

# Libraries to link with.
# For arcane reasons, statically linking with the pthread
# library on Linux requires "--whole-archive".
if(X86_64)
    target_link_libraries(
        shasta2Executable
        atomic boost_system boost_program_options boost_chrono boost_serialization  spoa png z
        lapack blas gfortran quadmath
        -Wl,--whole-archive -lpthread -Wl,--no-whole-archive)
else(X86_64)
    target_link_libraries(
        shasta2Executable
        atomic boost_system boost_program_options boost_chrono boost_serialization spoa png z
        lapack blas gfortran
        -Wl,--whole-archive -lpthread -Wl,--no-whole-archive)
endif(X86_64) 
  
# The static executable goes to the bin directory.
install(TARGETS shasta2Executable DESTINATION shasta-install/bin)


