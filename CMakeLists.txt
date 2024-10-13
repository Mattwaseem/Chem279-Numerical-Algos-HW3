cmake_minimum_required(VERSION 3.10)
project(MoleculeSetup)

set(CMAKE_CXX_STANDARD 17)

# Include directories
include_directories(include)

# Explicitly include Armadillo
set(ARMADILLO_INCLUDE_DIR /usr/local/opt/armadillo/include)  # Update to your Armadillo installation path
set(ARMADILLO_LIBRARY /usr/local/opt/armadillo/lib/libarmadillo.dylib)  # Update to your Armadillo library path

include_directories(${ARMADILLO_INCLUDE_DIR})

# Add source files
add_executable(MoleculeSetup
    src/main.cpp
    src/molecule.cpp
    src/input_parser.cpp
    src/CartesianGaussian.cpp
    src/OverlapMatrix.cpp
)

# Link Armadillo explicitly
target_link_libraries(MoleculeSetup ${ARMADILLO_LIBRARY})

# Output to bin
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/bin)

# "run" target that runs the program and generates the output file
add_custom_target(run
    COMMAND ${CMAKE_SOURCE_DIR}/bin/MoleculeSetup > ${CMAKE_SOURCE_DIR}/calculated_output_p1/C2H2_output.txt
    DEPENDS MoleculeSetup
    COMMENT "Running MoleculeSetup and generating output in calculated_output_p1/C2H2_output.txt"
    VERBATIM
)

# Custom clean-all target to remove build files and output
add_custom_target(clean-all
    COMMAND rm -rf ${CMAKE_BINARY_DIR}/*
    COMMAND rm -rf ${CMAKE_SOURCE_DIR}/bin/*
    COMMAND rm -rf ${CMAKE_SOURCE_DIR}/calculated_output_p1/*
    COMMENT "Removing all build files, executables, and output files"
)
