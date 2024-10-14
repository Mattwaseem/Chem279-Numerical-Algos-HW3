cmake_minimum_required(VERSION 3.10)
project(MoleculeSetup)

set(CMAKE_CXX_STANDARD 17)

# Include directories
include_directories(include)

# Explicitly set the Armadillo include and library paths
set(ARMADILLO_INCLUDE_DIR /opt/homebrew/opt/armadillo/include)
set(ARMADILLO_LIBRARY /opt/homebrew/opt/armadillo/lib/libarmadillo.dylib)

# Add the Armadillo include directory
include_directories(${ARMADILLO_INCLUDE_DIR})

# Add source files
add_executable(MoleculeSetup
    src/main.cpp
    src/molecule.cpp
    src/input_parser.cpp
    src/CartesianGaussian.cpp
    src/OverlapMatrix.cpp
)

# Link Armadillo
target_link_libraries(MoleculeSetup ${ARMADILLO_LIBRARY})

# Output to bin
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/bin)

# Run the program with different input files
add_custom_target(run
    COMMAND ${CMAKE_SOURCE_DIR}/bin/MoleculeSetup ${CMAKE_SOURCE_DIR}/sample_input/C2H2.txt > ${CMAKE_SOURCE_DIR}/calculated_output_p1/C2H2_output.txt
    COMMAND ${CMAKE_SOURCE_DIR}/bin/MoleculeSetup ${CMAKE_SOURCE_DIR}/sample_input/C2H2.txt > ${CMAKE_SOURCE_DIR}/calculated_output_p1/C2H4_output.txt
    COMMAND ${CMAKE_SOURCE_DIR}/bin/MoleculeSetup ${CMAKE_SOURCE_DIR}/sample_input/H2.txt > ${CMAKE_SOURCE_DIR}/calculated_output_p1/H2_output.txt
    DEPENDS MoleculeSetup
    COMMENT "Running MoleculeSetup with C2H2.txt and generating output"
    COMMENT "Running MoleculeSetup with H2.txt and generating output"
    COMMENT "Running MoleculeSetup with C2H4.txt and generating output"
    VERBATIM
)



# Custom clean-all target to remove build files and output
add_custom_target(clean-all
    COMMAND rm -rf ${CMAKE_BINARY_DIR}/*
    COMMAND rm -rf ${CMAKE_SOURCE_DIR}/bin/*
    COMMAND rm -rf ${CMAKE_SOURCE_DIR}/calculated_output_p1/*
    COMMENT "Removing all build files, executables, and output files"
)
