# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list

# Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /usr/bin/ccmake

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /scratch/data/anna/clonesim

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /scratch/data/anna/clonesim/build

# Include any dependencies generated for this target.
include CMakeFiles/simulate.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/simulate.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/simulate.dir/flags.make

CMakeFiles/simulate.dir/src/arg_parser.cpp.o: CMakeFiles/simulate.dir/flags.make
CMakeFiles/simulate.dir/src/arg_parser.cpp.o: ../src/arg_parser.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /scratch/data/anna/clonesim/build/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/simulate.dir/src/arg_parser.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/simulate.dir/src/arg_parser.cpp.o -c /scratch/data/anna/clonesim/src/arg_parser.cpp

CMakeFiles/simulate.dir/src/arg_parser.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/simulate.dir/src/arg_parser.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /scratch/data/anna/clonesim/src/arg_parser.cpp > CMakeFiles/simulate.dir/src/arg_parser.cpp.i

CMakeFiles/simulate.dir/src/arg_parser.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/simulate.dir/src/arg_parser.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /scratch/data/anna/clonesim/src/arg_parser.cpp -o CMakeFiles/simulate.dir/src/arg_parser.cpp.s

CMakeFiles/simulate.dir/src/arg_parser.cpp.o.requires:
.PHONY : CMakeFiles/simulate.dir/src/arg_parser.cpp.o.requires

CMakeFiles/simulate.dir/src/arg_parser.cpp.o.provides: CMakeFiles/simulate.dir/src/arg_parser.cpp.o.requires
	$(MAKE) -f CMakeFiles/simulate.dir/build.make CMakeFiles/simulate.dir/src/arg_parser.cpp.o.provides.build
.PHONY : CMakeFiles/simulate.dir/src/arg_parser.cpp.o.provides

CMakeFiles/simulate.dir/src/arg_parser.cpp.o.provides.build: CMakeFiles/simulate.dir/src/arg_parser.cpp.o

CMakeFiles/simulate.dir/src/simulatemain.cpp.o: CMakeFiles/simulate.dir/flags.make
CMakeFiles/simulate.dir/src/simulatemain.cpp.o: ../src/simulatemain.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /scratch/data/anna/clonesim/build/CMakeFiles $(CMAKE_PROGRESS_2)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/simulate.dir/src/simulatemain.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/simulate.dir/src/simulatemain.cpp.o -c /scratch/data/anna/clonesim/src/simulatemain.cpp

CMakeFiles/simulate.dir/src/simulatemain.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/simulate.dir/src/simulatemain.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /scratch/data/anna/clonesim/src/simulatemain.cpp > CMakeFiles/simulate.dir/src/simulatemain.cpp.i

CMakeFiles/simulate.dir/src/simulatemain.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/simulate.dir/src/simulatemain.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /scratch/data/anna/clonesim/src/simulatemain.cpp -o CMakeFiles/simulate.dir/src/simulatemain.cpp.s

CMakeFiles/simulate.dir/src/simulatemain.cpp.o.requires:
.PHONY : CMakeFiles/simulate.dir/src/simulatemain.cpp.o.requires

CMakeFiles/simulate.dir/src/simulatemain.cpp.o.provides: CMakeFiles/simulate.dir/src/simulatemain.cpp.o.requires
	$(MAKE) -f CMakeFiles/simulate.dir/build.make CMakeFiles/simulate.dir/src/simulatemain.cpp.o.provides.build
.PHONY : CMakeFiles/simulate.dir/src/simulatemain.cpp.o.provides

CMakeFiles/simulate.dir/src/simulatemain.cpp.o.provides.build: CMakeFiles/simulate.dir/src/simulatemain.cpp.o

CMakeFiles/simulate.dir/src/basetree.cpp.o: CMakeFiles/simulate.dir/flags.make
CMakeFiles/simulate.dir/src/basetree.cpp.o: ../src/basetree.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /scratch/data/anna/clonesim/build/CMakeFiles $(CMAKE_PROGRESS_3)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/simulate.dir/src/basetree.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/simulate.dir/src/basetree.cpp.o -c /scratch/data/anna/clonesim/src/basetree.cpp

CMakeFiles/simulate.dir/src/basetree.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/simulate.dir/src/basetree.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /scratch/data/anna/clonesim/src/basetree.cpp > CMakeFiles/simulate.dir/src/basetree.cpp.i

CMakeFiles/simulate.dir/src/basetree.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/simulate.dir/src/basetree.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /scratch/data/anna/clonesim/src/basetree.cpp -o CMakeFiles/simulate.dir/src/basetree.cpp.s

CMakeFiles/simulate.dir/src/basetree.cpp.o.requires:
.PHONY : CMakeFiles/simulate.dir/src/basetree.cpp.o.requires

CMakeFiles/simulate.dir/src/basetree.cpp.o.provides: CMakeFiles/simulate.dir/src/basetree.cpp.o.requires
	$(MAKE) -f CMakeFiles/simulate.dir/build.make CMakeFiles/simulate.dir/src/basetree.cpp.o.provides.build
.PHONY : CMakeFiles/simulate.dir/src/basetree.cpp.o.provides

CMakeFiles/simulate.dir/src/basetree.cpp.o.provides.build: CMakeFiles/simulate.dir/src/basetree.cpp.o

CMakeFiles/simulate.dir/src/cnatree.cpp.o: CMakeFiles/simulate.dir/flags.make
CMakeFiles/simulate.dir/src/cnatree.cpp.o: ../src/cnatree.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /scratch/data/anna/clonesim/build/CMakeFiles $(CMAKE_PROGRESS_4)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/simulate.dir/src/cnatree.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/simulate.dir/src/cnatree.cpp.o -c /scratch/data/anna/clonesim/src/cnatree.cpp

CMakeFiles/simulate.dir/src/cnatree.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/simulate.dir/src/cnatree.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /scratch/data/anna/clonesim/src/cnatree.cpp > CMakeFiles/simulate.dir/src/cnatree.cpp.i

CMakeFiles/simulate.dir/src/cnatree.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/simulate.dir/src/cnatree.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /scratch/data/anna/clonesim/src/cnatree.cpp -o CMakeFiles/simulate.dir/src/cnatree.cpp.s

CMakeFiles/simulate.dir/src/cnatree.cpp.o.requires:
.PHONY : CMakeFiles/simulate.dir/src/cnatree.cpp.o.requires

CMakeFiles/simulate.dir/src/cnatree.cpp.o.provides: CMakeFiles/simulate.dir/src/cnatree.cpp.o.requires
	$(MAKE) -f CMakeFiles/simulate.dir/build.make CMakeFiles/simulate.dir/src/cnatree.cpp.o.provides.build
.PHONY : CMakeFiles/simulate.dir/src/cnatree.cpp.o.provides

CMakeFiles/simulate.dir/src/cnatree.cpp.o.provides.build: CMakeFiles/simulate.dir/src/cnatree.cpp.o

CMakeFiles/simulate.dir/src/genotypetree.cpp.o: CMakeFiles/simulate.dir/flags.make
CMakeFiles/simulate.dir/src/genotypetree.cpp.o: ../src/genotypetree.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /scratch/data/anna/clonesim/build/CMakeFiles $(CMAKE_PROGRESS_5)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/simulate.dir/src/genotypetree.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/simulate.dir/src/genotypetree.cpp.o -c /scratch/data/anna/clonesim/src/genotypetree.cpp

CMakeFiles/simulate.dir/src/genotypetree.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/simulate.dir/src/genotypetree.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /scratch/data/anna/clonesim/src/genotypetree.cpp > CMakeFiles/simulate.dir/src/genotypetree.cpp.i

CMakeFiles/simulate.dir/src/genotypetree.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/simulate.dir/src/genotypetree.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /scratch/data/anna/clonesim/src/genotypetree.cpp -o CMakeFiles/simulate.dir/src/genotypetree.cpp.s

CMakeFiles/simulate.dir/src/genotypetree.cpp.o.requires:
.PHONY : CMakeFiles/simulate.dir/src/genotypetree.cpp.o.requires

CMakeFiles/simulate.dir/src/genotypetree.cpp.o.provides: CMakeFiles/simulate.dir/src/genotypetree.cpp.o.requires
	$(MAKE) -f CMakeFiles/simulate.dir/build.make CMakeFiles/simulate.dir/src/genotypetree.cpp.o.provides.build
.PHONY : CMakeFiles/simulate.dir/src/genotypetree.cpp.o.provides

CMakeFiles/simulate.dir/src/genotypetree.cpp.o.provides.build: CMakeFiles/simulate.dir/src/genotypetree.cpp.o

CMakeFiles/simulate.dir/src/genotypegraph.cpp.o: CMakeFiles/simulate.dir/flags.make
CMakeFiles/simulate.dir/src/genotypegraph.cpp.o: ../src/genotypegraph.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /scratch/data/anna/clonesim/build/CMakeFiles $(CMAKE_PROGRESS_6)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/simulate.dir/src/genotypegraph.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/simulate.dir/src/genotypegraph.cpp.o -c /scratch/data/anna/clonesim/src/genotypegraph.cpp

CMakeFiles/simulate.dir/src/genotypegraph.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/simulate.dir/src/genotypegraph.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /scratch/data/anna/clonesim/src/genotypegraph.cpp > CMakeFiles/simulate.dir/src/genotypegraph.cpp.i

CMakeFiles/simulate.dir/src/genotypegraph.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/simulate.dir/src/genotypegraph.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /scratch/data/anna/clonesim/src/genotypegraph.cpp -o CMakeFiles/simulate.dir/src/genotypegraph.cpp.s

CMakeFiles/simulate.dir/src/genotypegraph.cpp.o.requires:
.PHONY : CMakeFiles/simulate.dir/src/genotypegraph.cpp.o.requires

CMakeFiles/simulate.dir/src/genotypegraph.cpp.o.provides: CMakeFiles/simulate.dir/src/genotypegraph.cpp.o.requires
	$(MAKE) -f CMakeFiles/simulate.dir/build.make CMakeFiles/simulate.dir/src/genotypegraph.cpp.o.provides.build
.PHONY : CMakeFiles/simulate.dir/src/genotypegraph.cpp.o.provides

CMakeFiles/simulate.dir/src/genotypegraph.cpp.o.provides.build: CMakeFiles/simulate.dir/src/genotypegraph.cpp.o

CMakeFiles/simulate.dir/src/cnagraph.cpp.o: CMakeFiles/simulate.dir/flags.make
CMakeFiles/simulate.dir/src/cnagraph.cpp.o: ../src/cnagraph.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /scratch/data/anna/clonesim/build/CMakeFiles $(CMAKE_PROGRESS_7)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/simulate.dir/src/cnagraph.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/simulate.dir/src/cnagraph.cpp.o -c /scratch/data/anna/clonesim/src/cnagraph.cpp

CMakeFiles/simulate.dir/src/cnagraph.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/simulate.dir/src/cnagraph.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /scratch/data/anna/clonesim/src/cnagraph.cpp > CMakeFiles/simulate.dir/src/cnagraph.cpp.i

CMakeFiles/simulate.dir/src/cnagraph.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/simulate.dir/src/cnagraph.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /scratch/data/anna/clonesim/src/cnagraph.cpp -o CMakeFiles/simulate.dir/src/cnagraph.cpp.s

CMakeFiles/simulate.dir/src/cnagraph.cpp.o.requires:
.PHONY : CMakeFiles/simulate.dir/src/cnagraph.cpp.o.requires

CMakeFiles/simulate.dir/src/cnagraph.cpp.o.provides: CMakeFiles/simulate.dir/src/cnagraph.cpp.o.requires
	$(MAKE) -f CMakeFiles/simulate.dir/build.make CMakeFiles/simulate.dir/src/cnagraph.cpp.o.provides.build
.PHONY : CMakeFiles/simulate.dir/src/cnagraph.cpp.o.provides

CMakeFiles/simulate.dir/src/cnagraph.cpp.o.provides.build: CMakeFiles/simulate.dir/src/cnagraph.cpp.o

CMakeFiles/simulate.dir/src/utils.cpp.o: CMakeFiles/simulate.dir/flags.make
CMakeFiles/simulate.dir/src/utils.cpp.o: ../src/utils.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /scratch/data/anna/clonesim/build/CMakeFiles $(CMAKE_PROGRESS_8)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/simulate.dir/src/utils.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/simulate.dir/src/utils.cpp.o -c /scratch/data/anna/clonesim/src/utils.cpp

CMakeFiles/simulate.dir/src/utils.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/simulate.dir/src/utils.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /scratch/data/anna/clonesim/src/utils.cpp > CMakeFiles/simulate.dir/src/utils.cpp.i

CMakeFiles/simulate.dir/src/utils.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/simulate.dir/src/utils.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /scratch/data/anna/clonesim/src/utils.cpp -o CMakeFiles/simulate.dir/src/utils.cpp.s

CMakeFiles/simulate.dir/src/utils.cpp.o.requires:
.PHONY : CMakeFiles/simulate.dir/src/utils.cpp.o.requires

CMakeFiles/simulate.dir/src/utils.cpp.o.provides: CMakeFiles/simulate.dir/src/utils.cpp.o.requires
	$(MAKE) -f CMakeFiles/simulate.dir/build.make CMakeFiles/simulate.dir/src/utils.cpp.o.provides.build
.PHONY : CMakeFiles/simulate.dir/src/utils.cpp.o.provides

CMakeFiles/simulate.dir/src/utils.cpp.o.provides.build: CMakeFiles/simulate.dir/src/utils.cpp.o

CMakeFiles/simulate.dir/src/basematrix.cpp.o: CMakeFiles/simulate.dir/flags.make
CMakeFiles/simulate.dir/src/basematrix.cpp.o: ../src/basematrix.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /scratch/data/anna/clonesim/build/CMakeFiles $(CMAKE_PROGRESS_9)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/simulate.dir/src/basematrix.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/simulate.dir/src/basematrix.cpp.o -c /scratch/data/anna/clonesim/src/basematrix.cpp

CMakeFiles/simulate.dir/src/basematrix.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/simulate.dir/src/basematrix.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /scratch/data/anna/clonesim/src/basematrix.cpp > CMakeFiles/simulate.dir/src/basematrix.cpp.i

CMakeFiles/simulate.dir/src/basematrix.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/simulate.dir/src/basematrix.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /scratch/data/anna/clonesim/src/basematrix.cpp -o CMakeFiles/simulate.dir/src/basematrix.cpp.s

CMakeFiles/simulate.dir/src/basematrix.cpp.o.requires:
.PHONY : CMakeFiles/simulate.dir/src/basematrix.cpp.o.requires

CMakeFiles/simulate.dir/src/basematrix.cpp.o.provides: CMakeFiles/simulate.dir/src/basematrix.cpp.o.requires
	$(MAKE) -f CMakeFiles/simulate.dir/build.make CMakeFiles/simulate.dir/src/basematrix.cpp.o.provides.build
.PHONY : CMakeFiles/simulate.dir/src/basematrix.cpp.o.provides

CMakeFiles/simulate.dir/src/basematrix.cpp.o.provides.build: CMakeFiles/simulate.dir/src/basematrix.cpp.o

CMakeFiles/simulate.dir/src/phylogeny.cpp.o: CMakeFiles/simulate.dir/flags.make
CMakeFiles/simulate.dir/src/phylogeny.cpp.o: ../src/phylogeny.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /scratch/data/anna/clonesim/build/CMakeFiles $(CMAKE_PROGRESS_10)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/simulate.dir/src/phylogeny.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/simulate.dir/src/phylogeny.cpp.o -c /scratch/data/anna/clonesim/src/phylogeny.cpp

CMakeFiles/simulate.dir/src/phylogeny.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/simulate.dir/src/phylogeny.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /scratch/data/anna/clonesim/src/phylogeny.cpp > CMakeFiles/simulate.dir/src/phylogeny.cpp.i

CMakeFiles/simulate.dir/src/phylogeny.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/simulate.dir/src/phylogeny.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /scratch/data/anna/clonesim/src/phylogeny.cpp -o CMakeFiles/simulate.dir/src/phylogeny.cpp.s

CMakeFiles/simulate.dir/src/phylogeny.cpp.o.requires:
.PHONY : CMakeFiles/simulate.dir/src/phylogeny.cpp.o.requires

CMakeFiles/simulate.dir/src/phylogeny.cpp.o.provides: CMakeFiles/simulate.dir/src/phylogeny.cpp.o.requires
	$(MAKE) -f CMakeFiles/simulate.dir/build.make CMakeFiles/simulate.dir/src/phylogeny.cpp.o.provides.build
.PHONY : CMakeFiles/simulate.dir/src/phylogeny.cpp.o.provides

CMakeFiles/simulate.dir/src/phylogeny.cpp.o.provides.build: CMakeFiles/simulate.dir/src/phylogeny.cpp.o

# Object files for target simulate
simulate_OBJECTS = \
"CMakeFiles/simulate.dir/src/arg_parser.cpp.o" \
"CMakeFiles/simulate.dir/src/simulatemain.cpp.o" \
"CMakeFiles/simulate.dir/src/basetree.cpp.o" \
"CMakeFiles/simulate.dir/src/cnatree.cpp.o" \
"CMakeFiles/simulate.dir/src/genotypetree.cpp.o" \
"CMakeFiles/simulate.dir/src/genotypegraph.cpp.o" \
"CMakeFiles/simulate.dir/src/cnagraph.cpp.o" \
"CMakeFiles/simulate.dir/src/utils.cpp.o" \
"CMakeFiles/simulate.dir/src/basematrix.cpp.o" \
"CMakeFiles/simulate.dir/src/phylogeny.cpp.o"

# External object files for target simulate
simulate_EXTERNAL_OBJECTS =

simulate: CMakeFiles/simulate.dir/src/arg_parser.cpp.o
simulate: CMakeFiles/simulate.dir/src/simulatemain.cpp.o
simulate: CMakeFiles/simulate.dir/src/basetree.cpp.o
simulate: CMakeFiles/simulate.dir/src/cnatree.cpp.o
simulate: CMakeFiles/simulate.dir/src/genotypetree.cpp.o
simulate: CMakeFiles/simulate.dir/src/genotypegraph.cpp.o
simulate: CMakeFiles/simulate.dir/src/cnagraph.cpp.o
simulate: CMakeFiles/simulate.dir/src/utils.cpp.o
simulate: CMakeFiles/simulate.dir/src/basematrix.cpp.o
simulate: CMakeFiles/simulate.dir/src/phylogeny.cpp.o
simulate: CMakeFiles/simulate.dir/build.make
simulate: /usr/lib64/libboost_thread-mt.so
simulate: /usr/lib64/libboost_system-mt.so
simulate: /usr/lib64/libboost_filesystem-mt.so
simulate: CMakeFiles/simulate.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable simulate"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/simulate.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/simulate.dir/build: simulate
.PHONY : CMakeFiles/simulate.dir/build

CMakeFiles/simulate.dir/requires: CMakeFiles/simulate.dir/src/arg_parser.cpp.o.requires
CMakeFiles/simulate.dir/requires: CMakeFiles/simulate.dir/src/simulatemain.cpp.o.requires
CMakeFiles/simulate.dir/requires: CMakeFiles/simulate.dir/src/basetree.cpp.o.requires
CMakeFiles/simulate.dir/requires: CMakeFiles/simulate.dir/src/cnatree.cpp.o.requires
CMakeFiles/simulate.dir/requires: CMakeFiles/simulate.dir/src/genotypetree.cpp.o.requires
CMakeFiles/simulate.dir/requires: CMakeFiles/simulate.dir/src/genotypegraph.cpp.o.requires
CMakeFiles/simulate.dir/requires: CMakeFiles/simulate.dir/src/cnagraph.cpp.o.requires
CMakeFiles/simulate.dir/requires: CMakeFiles/simulate.dir/src/utils.cpp.o.requires
CMakeFiles/simulate.dir/requires: CMakeFiles/simulate.dir/src/basematrix.cpp.o.requires
CMakeFiles/simulate.dir/requires: CMakeFiles/simulate.dir/src/phylogeny.cpp.o.requires
.PHONY : CMakeFiles/simulate.dir/requires

CMakeFiles/simulate.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/simulate.dir/cmake_clean.cmake
.PHONY : CMakeFiles/simulate.dir/clean

CMakeFiles/simulate.dir/depend:
	cd /scratch/data/anna/clonesim/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /scratch/data/anna/clonesim /scratch/data/anna/clonesim /scratch/data/anna/clonesim/build /scratch/data/anna/clonesim/build /scratch/data/anna/clonesim/build/CMakeFiles/simulate.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/simulate.dir/depend

