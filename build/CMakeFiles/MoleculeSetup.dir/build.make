# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.30

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /opt/homebrew/Cellar/cmake/3.30.3/bin/cmake

# The command to remove a file.
RM = /opt/homebrew/Cellar/cmake/3.30.3/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/nematsfolder/chem279/homework3

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/nematsfolder/chem279/homework3/build

# Include any dependencies generated for this target.
include CMakeFiles/MoleculeSetup.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/MoleculeSetup.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/MoleculeSetup.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/MoleculeSetup.dir/flags.make

CMakeFiles/MoleculeSetup.dir/src/main.cpp.o: CMakeFiles/MoleculeSetup.dir/flags.make
CMakeFiles/MoleculeSetup.dir/src/main.cpp.o: /Users/nematsfolder/chem279/homework3/src/main.cpp
CMakeFiles/MoleculeSetup.dir/src/main.cpp.o: CMakeFiles/MoleculeSetup.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/nematsfolder/chem279/homework3/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/MoleculeSetup.dir/src/main.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/MoleculeSetup.dir/src/main.cpp.o -MF CMakeFiles/MoleculeSetup.dir/src/main.cpp.o.d -o CMakeFiles/MoleculeSetup.dir/src/main.cpp.o -c /Users/nematsfolder/chem279/homework3/src/main.cpp

CMakeFiles/MoleculeSetup.dir/src/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/MoleculeSetup.dir/src/main.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/nematsfolder/chem279/homework3/src/main.cpp > CMakeFiles/MoleculeSetup.dir/src/main.cpp.i

CMakeFiles/MoleculeSetup.dir/src/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/MoleculeSetup.dir/src/main.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/nematsfolder/chem279/homework3/src/main.cpp -o CMakeFiles/MoleculeSetup.dir/src/main.cpp.s

CMakeFiles/MoleculeSetup.dir/src/molecule.cpp.o: CMakeFiles/MoleculeSetup.dir/flags.make
CMakeFiles/MoleculeSetup.dir/src/molecule.cpp.o: /Users/nematsfolder/chem279/homework3/src/molecule.cpp
CMakeFiles/MoleculeSetup.dir/src/molecule.cpp.o: CMakeFiles/MoleculeSetup.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/nematsfolder/chem279/homework3/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/MoleculeSetup.dir/src/molecule.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/MoleculeSetup.dir/src/molecule.cpp.o -MF CMakeFiles/MoleculeSetup.dir/src/molecule.cpp.o.d -o CMakeFiles/MoleculeSetup.dir/src/molecule.cpp.o -c /Users/nematsfolder/chem279/homework3/src/molecule.cpp

CMakeFiles/MoleculeSetup.dir/src/molecule.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/MoleculeSetup.dir/src/molecule.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/nematsfolder/chem279/homework3/src/molecule.cpp > CMakeFiles/MoleculeSetup.dir/src/molecule.cpp.i

CMakeFiles/MoleculeSetup.dir/src/molecule.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/MoleculeSetup.dir/src/molecule.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/nematsfolder/chem279/homework3/src/molecule.cpp -o CMakeFiles/MoleculeSetup.dir/src/molecule.cpp.s

CMakeFiles/MoleculeSetup.dir/src/input_parser.cpp.o: CMakeFiles/MoleculeSetup.dir/flags.make
CMakeFiles/MoleculeSetup.dir/src/input_parser.cpp.o: /Users/nematsfolder/chem279/homework3/src/input_parser.cpp
CMakeFiles/MoleculeSetup.dir/src/input_parser.cpp.o: CMakeFiles/MoleculeSetup.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/nematsfolder/chem279/homework3/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/MoleculeSetup.dir/src/input_parser.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/MoleculeSetup.dir/src/input_parser.cpp.o -MF CMakeFiles/MoleculeSetup.dir/src/input_parser.cpp.o.d -o CMakeFiles/MoleculeSetup.dir/src/input_parser.cpp.o -c /Users/nematsfolder/chem279/homework3/src/input_parser.cpp

CMakeFiles/MoleculeSetup.dir/src/input_parser.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/MoleculeSetup.dir/src/input_parser.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/nematsfolder/chem279/homework3/src/input_parser.cpp > CMakeFiles/MoleculeSetup.dir/src/input_parser.cpp.i

CMakeFiles/MoleculeSetup.dir/src/input_parser.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/MoleculeSetup.dir/src/input_parser.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/nematsfolder/chem279/homework3/src/input_parser.cpp -o CMakeFiles/MoleculeSetup.dir/src/input_parser.cpp.s

CMakeFiles/MoleculeSetup.dir/src/CartesianGaussian.cpp.o: CMakeFiles/MoleculeSetup.dir/flags.make
CMakeFiles/MoleculeSetup.dir/src/CartesianGaussian.cpp.o: /Users/nematsfolder/chem279/homework3/src/CartesianGaussian.cpp
CMakeFiles/MoleculeSetup.dir/src/CartesianGaussian.cpp.o: CMakeFiles/MoleculeSetup.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/nematsfolder/chem279/homework3/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/MoleculeSetup.dir/src/CartesianGaussian.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/MoleculeSetup.dir/src/CartesianGaussian.cpp.o -MF CMakeFiles/MoleculeSetup.dir/src/CartesianGaussian.cpp.o.d -o CMakeFiles/MoleculeSetup.dir/src/CartesianGaussian.cpp.o -c /Users/nematsfolder/chem279/homework3/src/CartesianGaussian.cpp

CMakeFiles/MoleculeSetup.dir/src/CartesianGaussian.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/MoleculeSetup.dir/src/CartesianGaussian.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/nematsfolder/chem279/homework3/src/CartesianGaussian.cpp > CMakeFiles/MoleculeSetup.dir/src/CartesianGaussian.cpp.i

CMakeFiles/MoleculeSetup.dir/src/CartesianGaussian.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/MoleculeSetup.dir/src/CartesianGaussian.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/nematsfolder/chem279/homework3/src/CartesianGaussian.cpp -o CMakeFiles/MoleculeSetup.dir/src/CartesianGaussian.cpp.s

CMakeFiles/MoleculeSetup.dir/src/OverlapMatrix.cpp.o: CMakeFiles/MoleculeSetup.dir/flags.make
CMakeFiles/MoleculeSetup.dir/src/OverlapMatrix.cpp.o: /Users/nematsfolder/chem279/homework3/src/OverlapMatrix.cpp
CMakeFiles/MoleculeSetup.dir/src/OverlapMatrix.cpp.o: CMakeFiles/MoleculeSetup.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/nematsfolder/chem279/homework3/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMakeFiles/MoleculeSetup.dir/src/OverlapMatrix.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/MoleculeSetup.dir/src/OverlapMatrix.cpp.o -MF CMakeFiles/MoleculeSetup.dir/src/OverlapMatrix.cpp.o.d -o CMakeFiles/MoleculeSetup.dir/src/OverlapMatrix.cpp.o -c /Users/nematsfolder/chem279/homework3/src/OverlapMatrix.cpp

CMakeFiles/MoleculeSetup.dir/src/OverlapMatrix.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/MoleculeSetup.dir/src/OverlapMatrix.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/nematsfolder/chem279/homework3/src/OverlapMatrix.cpp > CMakeFiles/MoleculeSetup.dir/src/OverlapMatrix.cpp.i

CMakeFiles/MoleculeSetup.dir/src/OverlapMatrix.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/MoleculeSetup.dir/src/OverlapMatrix.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/nematsfolder/chem279/homework3/src/OverlapMatrix.cpp -o CMakeFiles/MoleculeSetup.dir/src/OverlapMatrix.cpp.s

# Object files for target MoleculeSetup
MoleculeSetup_OBJECTS = \
"CMakeFiles/MoleculeSetup.dir/src/main.cpp.o" \
"CMakeFiles/MoleculeSetup.dir/src/molecule.cpp.o" \
"CMakeFiles/MoleculeSetup.dir/src/input_parser.cpp.o" \
"CMakeFiles/MoleculeSetup.dir/src/CartesianGaussian.cpp.o" \
"CMakeFiles/MoleculeSetup.dir/src/OverlapMatrix.cpp.o"

# External object files for target MoleculeSetup
MoleculeSetup_EXTERNAL_OBJECTS =

/Users/nematsfolder/chem279/homework3/bin/MoleculeSetup: CMakeFiles/MoleculeSetup.dir/src/main.cpp.o
/Users/nematsfolder/chem279/homework3/bin/MoleculeSetup: CMakeFiles/MoleculeSetup.dir/src/molecule.cpp.o
/Users/nematsfolder/chem279/homework3/bin/MoleculeSetup: CMakeFiles/MoleculeSetup.dir/src/input_parser.cpp.o
/Users/nematsfolder/chem279/homework3/bin/MoleculeSetup: CMakeFiles/MoleculeSetup.dir/src/CartesianGaussian.cpp.o
/Users/nematsfolder/chem279/homework3/bin/MoleculeSetup: CMakeFiles/MoleculeSetup.dir/src/OverlapMatrix.cpp.o
/Users/nematsfolder/chem279/homework3/bin/MoleculeSetup: CMakeFiles/MoleculeSetup.dir/build.make
/Users/nematsfolder/chem279/homework3/bin/MoleculeSetup: CMakeFiles/MoleculeSetup.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=/Users/nematsfolder/chem279/homework3/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Linking CXX executable /Users/nematsfolder/chem279/homework3/bin/MoleculeSetup"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/MoleculeSetup.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/MoleculeSetup.dir/build: /Users/nematsfolder/chem279/homework3/bin/MoleculeSetup
.PHONY : CMakeFiles/MoleculeSetup.dir/build

CMakeFiles/MoleculeSetup.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/MoleculeSetup.dir/cmake_clean.cmake
.PHONY : CMakeFiles/MoleculeSetup.dir/clean

CMakeFiles/MoleculeSetup.dir/depend:
	cd /Users/nematsfolder/chem279/homework3/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/nematsfolder/chem279/homework3 /Users/nematsfolder/chem279/homework3 /Users/nematsfolder/chem279/homework3/build /Users/nematsfolder/chem279/homework3/build /Users/nematsfolder/chem279/homework3/build/CMakeFiles/MoleculeSetup.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : CMakeFiles/MoleculeSetup.dir/depend

