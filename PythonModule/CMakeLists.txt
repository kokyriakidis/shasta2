cmake_minimum_required(VERSION 3.16)

project(shasta2PythonModule)

# C++ dialect.
add_definitions(-std=c++23)

# Compilation warnings.
add_definitions(-Wall -Wconversion -Wno-unused-result -Wno-trigraphs -Wno-psabi -Wunused-parameter -Wno-unused-function)

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

# This makes sure the code in PythonModule.cpp is compiled.
add_definitions(-DSHASTA_PYTHON_API)

# Sources files.
file(GLOB SOURCES ../src/*.cpp)

# Include directory.
include_directories(../src ~/abpoa/abPOA/include)

# Define our library.
add_library(shasta2PythonModule SHARED ${SOURCES})

# Make sure the library is named shasta.so.
set_target_properties(shasta2PythonModule PROPERTIES OUTPUT_NAME "shasta2")
set_target_properties(shasta2PythonModule PROPERTIES PREFIX "")
set_target_properties(shasta2PythonModule PROPERTIES DEFINE_SYMBOL "")



# Include directories and libraries for Python3 and pybind11.
execute_process(COMMAND /usr/bin/python3-config --embed --libs OUTPUT_VARIABLE SHASTA2_PYTHON_LIBRARIES)
execute_process(COMMAND python3 -m pybind11 --includes OUTPUT_VARIABLE SHASTA2_PYTHON_INCLUDES)
add_definitions(${SHASTA2_PYTHON_INCLUDES})
string(STRIP ${SHASTA2_PYTHON_LIBRARIES} SHASTA2_PYTHON_LIBRARIES)
SET(CMAKE_LINKER_FLAGS  "${CMAKE_LINKER_FLAGS} ${SHASTA2_PYTHON_LIBRARIES}")



# Libraries to link with.
target_link_libraries(
    shasta2PythonModule 
    atomic png boost_program_options boost_serialization pthread abpoa z m spoa lapack blas ${SHASTA2_PYTHON_LIBRARIES})

# Install the shared library into the bin directory.
install(TARGETS shasta2PythonModule DESTINATION shasta2-install/bin)




