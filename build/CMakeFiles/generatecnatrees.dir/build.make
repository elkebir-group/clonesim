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
include CMakeFiles/generatecnatrees.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/generatecnatrees.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/generatecnatrees.dir/flags.make

CMakeFiles/generatecnatrees.dir/src/arg_parser.cpp.o: CMakeFiles/generatecnatrees.dir/flags.make
CMakeFiles/generatecnatrees.dir/src/arg_parser.cpp.o: ../src/arg_parser.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /scratch/data/anna/clonesim/build/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/generatecnatrees.dir/src/arg_parser.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/generatecnatrees.dir/src/arg_parser.cpp.o -c /scratch/data/anna/clonesim/src/arg_parser.cpp

CMakeFiles/generatecnatrees.dir/src/arg_parser.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/generatecnatrees.dir/src/arg_parser.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /scratch/data/anna/clonesim/src/arg_parser.cpp > CMakeFiles/generatecnatrees.dir/src/arg_parser.cpp.i

CMakeFiles/generatecnatrees.dir/src/arg_parser.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/generatecnatrees.dir/src/arg_parser.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /scratch/data/anna/clonesim/src/arg_parser.cpp -o CMakeFiles/generatecnatrees.dir/src/arg_parser.cpp.s

CMakeFiles/generatecnatrees.dir/src/arg_parser.cpp.o.requires:
.PHONY : CMakeFiles/generatecnatrees.dir/src/arg_parser.cpp.o.requires

CMakeFiles/generatecnatrees.dir/src/arg_parser.cpp.o.provides: CMakeFiles/generatecnatrees.dir/src/arg_parser.cpp.o.requires
	$(MAKE) -f CMakeFiles/generatecnatrees.dir/build.make CMakeFiles/generatecnatrees.dir/src/arg_parser.cpp.o.provides.build
.PHONY : CMakeFiles/generatecnatrees.dir/src/arg_parser.cpp.o.provides

CMakeFiles/generatecnatrees.dir/src/arg_parser.cpp.o.provides.build: CMakeFiles/generatecnatrees.dir/src/arg_parser.cpp.o

CMakeFiles/generatecnatrees.dir/src/generatecnatreesmain.cpp.o: CMakeFiles/generatecnatrees.dir/flags.make
CMakeFiles/generatecnatrees.dir/src/generatecnatreesmain.cpp.o: ../src/generatecnatreesmain.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /scratch/data/anna/clonesim/build/CMakeFiles $(CMAKE_PROGRESS_2)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/generatecnatrees.dir/src/generatecnatreesmain.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/generatecnatrees.dir/src/generatecnatreesmain.cpp.o -c /scratch/data/anna/clonesim/src/generatecnatreesmain.cpp

CMakeFiles/generatecnatrees.dir/src/generatecnatreesmain.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/generatecnatrees.dir/src/generatecnatreesmain.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /scratch/data/anna/clonesim/src/generatecnatreesmain.cpp > CMakeFiles/generatecnatrees.dir/src/generatecnatreesmain.cpp.i

CMakeFiles/generatecnatrees.dir/src/generatecnatreesmain.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/generatecnatrees.dir/src/generatecnatreesmain.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /scratch/data/anna/clonesim/src/generatecnatreesmain.cpp -o CMakeFiles/generatecnatrees.dir/src/generatecnatreesmain.cpp.s

CMakeFiles/generatecnatrees.dir/src/generatecnatreesmain.cpp.o.requires:
.PHONY : CMakeFiles/generatecnatrees.dir/src/generatecnatreesmain.cpp.o.requires

CMakeFiles/generatecnatrees.dir/src/generatecnatreesmain.cpp.o.provides: CMakeFiles/generatecnatrees.dir/src/generatecnatreesmain.cpp.o.requires
	$(MAKE) -f CMakeFiles/generatecnatrees.dir/build.make CMakeFiles/generatecnatrees.dir/src/generatecnatreesmain.cpp.o.provides.build
.PHONY : CMakeFiles/generatecnatrees.dir/src/generatecnatreesmain.cpp.o.provides

CMakeFiles/generatecnatrees.dir/src/generatecnatreesmain.cpp.o.provides.build: CMakeFiles/generatecnatrees.dir/src/generatecnatreesmain.cpp.o

CMakeFiles/generatecnatrees.dir/src/cnatree.cpp.o: CMakeFiles/generatecnatrees.dir/flags.make
CMakeFiles/generatecnatrees.dir/src/cnatree.cpp.o: ../src/cnatree.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /scratch/data/anna/clonesim/build/CMakeFiles $(CMAKE_PROGRESS_3)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/generatecnatrees.dir/src/cnatree.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/generatecnatrees.dir/src/cnatree.cpp.o -c /scratch/data/anna/clonesim/src/cnatree.cpp

CMakeFiles/generatecnatrees.dir/src/cnatree.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/generatecnatrees.dir/src/cnatree.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /scratch/data/anna/clonesim/src/cnatree.cpp > CMakeFiles/generatecnatrees.dir/src/cnatree.cpp.i

CMakeFiles/generatecnatrees.dir/src/cnatree.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/generatecnatrees.dir/src/cnatree.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /scratch/data/anna/clonesim/src/cnatree.cpp -o CMakeFiles/generatecnatrees.dir/src/cnatree.cpp.s

CMakeFiles/generatecnatrees.dir/src/cnatree.cpp.o.requires:
.PHONY : CMakeFiles/generatecnatrees.dir/src/cnatree.cpp.o.requires

CMakeFiles/generatecnatrees.dir/src/cnatree.cpp.o.provides: CMakeFiles/generatecnatrees.dir/src/cnatree.cpp.o.requires
	$(MAKE) -f CMakeFiles/generatecnatrees.dir/build.make CMakeFiles/generatecnatrees.dir/src/cnatree.cpp.o.provides.build
.PHONY : CMakeFiles/generatecnatrees.dir/src/cnatree.cpp.o.provides

CMakeFiles/generatecnatrees.dir/src/cnatree.cpp.o.provides.build: CMakeFiles/generatecnatrees.dir/src/cnatree.cpp.o

CMakeFiles/generatecnatrees.dir/src/cnagraph.cpp.o: CMakeFiles/generatecnatrees.dir/flags.make
CMakeFiles/generatecnatrees.dir/src/cnagraph.cpp.o: ../src/cnagraph.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /scratch/data/anna/clonesim/build/CMakeFiles $(CMAKE_PROGRESS_4)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/generatecnatrees.dir/src/cnagraph.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/generatecnatrees.dir/src/cnagraph.cpp.o -c /scratch/data/anna/clonesim/src/cnagraph.cpp

CMakeFiles/generatecnatrees.dir/src/cnagraph.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/generatecnatrees.dir/src/cnagraph.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /scratch/data/anna/clonesim/src/cnagraph.cpp > CMakeFiles/generatecnatrees.dir/src/cnagraph.cpp.i

CMakeFiles/generatecnatrees.dir/src/cnagraph.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/generatecnatrees.dir/src/cnagraph.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /scratch/data/anna/clonesim/src/cnagraph.cpp -o CMakeFiles/generatecnatrees.dir/src/cnagraph.cpp.s

CMakeFiles/generatecnatrees.dir/src/cnagraph.cpp.o.requires:
.PHONY : CMakeFiles/generatecnatrees.dir/src/cnagraph.cpp.o.requires

CMakeFiles/generatecnatrees.dir/src/cnagraph.cpp.o.provides: CMakeFiles/generatecnatrees.dir/src/cnagraph.cpp.o.requires
	$(MAKE) -f CMakeFiles/generatecnatrees.dir/build.make CMakeFiles/generatecnatrees.dir/src/cnagraph.cpp.o.provides.build
.PHONY : CMakeFiles/generatecnatrees.dir/src/cnagraph.cpp.o.provides

CMakeFiles/generatecnatrees.dir/src/cnagraph.cpp.o.provides.build: CMakeFiles/generatecnatrees.dir/src/cnagraph.cpp.o

CMakeFiles/generatecnatrees.dir/src/genotypetree.cpp.o: CMakeFiles/generatecnatrees.dir/flags.make
CMakeFiles/generatecnatrees.dir/src/genotypetree.cpp.o: ../src/genotypetree.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /scratch/data/anna/clonesim/build/CMakeFiles $(CMAKE_PROGRESS_5)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/generatecnatrees.dir/src/genotypetree.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/generatecnatrees.dir/src/genotypetree.cpp.o -c /scratch/data/anna/clonesim/src/genotypetree.cpp

CMakeFiles/generatecnatrees.dir/src/genotypetree.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/generatecnatrees.dir/src/genotypetree.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /scratch/data/anna/clonesim/src/genotypetree.cpp > CMakeFiles/generatecnatrees.dir/src/genotypetree.cpp.i

CMakeFiles/generatecnatrees.dir/src/genotypetree.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/generatecnatrees.dir/src/genotypetree.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /scratch/data/anna/clonesim/src/genotypetree.cpp -o CMakeFiles/generatecnatrees.dir/src/genotypetree.cpp.s

CMakeFiles/generatecnatrees.dir/src/genotypetree.cpp.o.requires:
.PHONY : CMakeFiles/generatecnatrees.dir/src/genotypetree.cpp.o.requires

CMakeFiles/generatecnatrees.dir/src/genotypetree.cpp.o.provides: CMakeFiles/generatecnatrees.dir/src/genotypetree.cpp.o.requires
	$(MAKE) -f CMakeFiles/generatecnatrees.dir/build.make CMakeFiles/generatecnatrees.dir/src/genotypetree.cpp.o.provides.build
.PHONY : CMakeFiles/generatecnatrees.dir/src/genotypetree.cpp.o.provides

CMakeFiles/generatecnatrees.dir/src/genotypetree.cpp.o.provides.build: CMakeFiles/generatecnatrees.dir/src/genotypetree.cpp.o

CMakeFiles/generatecnatrees.dir/src/utils.cpp.o: CMakeFiles/generatecnatrees.dir/flags.make
CMakeFiles/generatecnatrees.dir/src/utils.cpp.o: ../src/utils.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /scratch/data/anna/clonesim/build/CMakeFiles $(CMAKE_PROGRESS_6)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/generatecnatrees.dir/src/utils.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/generatecnatrees.dir/src/utils.cpp.o -c /scratch/data/anna/clonesim/src/utils.cpp

CMakeFiles/generatecnatrees.dir/src/utils.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/generatecnatrees.dir/src/utils.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /scratch/data/anna/clonesim/src/utils.cpp > CMakeFiles/generatecnatrees.dir/src/utils.cpp.i

CMakeFiles/generatecnatrees.dir/src/utils.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/generatecnatrees.dir/src/utils.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /scratch/data/anna/clonesim/src/utils.cpp -o CMakeFiles/generatecnatrees.dir/src/utils.cpp.s

CMakeFiles/generatecnatrees.dir/src/utils.cpp.o.requires:
.PHONY : CMakeFiles/generatecnatrees.dir/src/utils.cpp.o.requires

CMakeFiles/generatecnatrees.dir/src/utils.cpp.o.provides: CMakeFiles/generatecnatrees.dir/src/utils.cpp.o.requires
	$(MAKE) -f CMakeFiles/generatecnatrees.dir/build.make CMakeFiles/generatecnatrees.dir/src/utils.cpp.o.provides.build
.PHONY : CMakeFiles/generatecnatrees.dir/src/utils.cpp.o.provides

CMakeFiles/generatecnatrees.dir/src/utils.cpp.o.provides.build: CMakeFiles/generatecnatrees.dir/src/utils.cpp.o

CMakeFiles/generatecnatrees.dir/src/basematrix.cpp.o: CMakeFiles/generatecnatrees.dir/flags.make
CMakeFiles/generatecnatrees.dir/src/basematrix.cpp.o: ../src/basematrix.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /scratch/data/anna/clonesim/build/CMakeFiles $(CMAKE_PROGRESS_7)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/generatecnatrees.dir/src/basematrix.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/generatecnatrees.dir/src/basematrix.cpp.o -c /scratch/data/anna/clonesim/src/basematrix.cpp

CMakeFiles/generatecnatrees.dir/src/basematrix.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/generatecnatrees.dir/src/basematrix.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /scratch/data/anna/clonesim/src/basematrix.cpp > CMakeFiles/generatecnatrees.dir/src/basematrix.cpp.i

CMakeFiles/generatecnatrees.dir/src/basematrix.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/generatecnatrees.dir/src/basematrix.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /scratch/data/anna/clonesim/src/basematrix.cpp -o CMakeFiles/generatecnatrees.dir/src/basematrix.cpp.s

CMakeFiles/generatecnatrees.dir/src/basematrix.cpp.o.requires:
.PHONY : CMakeFiles/generatecnatrees.dir/src/basematrix.cpp.o.requires

CMakeFiles/generatecnatrees.dir/src/basematrix.cpp.o.provides: CMakeFiles/generatecnatrees.dir/src/basematrix.cpp.o.requires
	$(MAKE) -f CMakeFiles/generatecnatrees.dir/build.make CMakeFiles/generatecnatrees.dir/src/basematrix.cpp.o.provides.build
.PHONY : CMakeFiles/generatecnatrees.dir/src/basematrix.cpp.o.provides

CMakeFiles/generatecnatrees.dir/src/basematrix.cpp.o.provides.build: CMakeFiles/generatecnatrees.dir/src/basematrix.cpp.o

CMakeFiles/generatecnatrees.dir/src/basetree.cpp.o: CMakeFiles/generatecnatrees.dir/flags.make
CMakeFiles/generatecnatrees.dir/src/basetree.cpp.o: ../src/basetree.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /scratch/data/anna/clonesim/build/CMakeFiles $(CMAKE_PROGRESS_8)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/generatecnatrees.dir/src/basetree.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/generatecnatrees.dir/src/basetree.cpp.o -c /scratch/data/anna/clonesim/src/basetree.cpp

CMakeFiles/generatecnatrees.dir/src/basetree.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/generatecnatrees.dir/src/basetree.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /scratch/data/anna/clonesim/src/basetree.cpp > CMakeFiles/generatecnatrees.dir/src/basetree.cpp.i

CMakeFiles/generatecnatrees.dir/src/basetree.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/generatecnatrees.dir/src/basetree.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /scratch/data/anna/clonesim/src/basetree.cpp -o CMakeFiles/generatecnatrees.dir/src/basetree.cpp.s

CMakeFiles/generatecnatrees.dir/src/basetree.cpp.o.requires:
.PHONY : CMakeFiles/generatecnatrees.dir/src/basetree.cpp.o.requires

CMakeFiles/generatecnatrees.dir/src/basetree.cpp.o.provides: CMakeFiles/generatecnatrees.dir/src/basetree.cpp.o.requires
	$(MAKE) -f CMakeFiles/generatecnatrees.dir/build.make CMakeFiles/generatecnatrees.dir/src/basetree.cpp.o.provides.build
.PHONY : CMakeFiles/generatecnatrees.dir/src/basetree.cpp.o.provides

CMakeFiles/generatecnatrees.dir/src/basetree.cpp.o.provides.build: CMakeFiles/generatecnatrees.dir/src/basetree.cpp.o

# Object files for target generatecnatrees
generatecnatrees_OBJECTS = \
"CMakeFiles/generatecnatrees.dir/src/arg_parser.cpp.o" \
"CMakeFiles/generatecnatrees.dir/src/generatecnatreesmain.cpp.o" \
"CMakeFiles/generatecnatrees.dir/src/cnatree.cpp.o" \
"CMakeFiles/generatecnatrees.dir/src/cnagraph.cpp.o" \
"CMakeFiles/generatecnatrees.dir/src/genotypetree.cpp.o" \
"CMakeFiles/generatecnatrees.dir/src/utils.cpp.o" \
"CMakeFiles/generatecnatrees.dir/src/basematrix.cpp.o" \
"CMakeFiles/generatecnatrees.dir/src/basetree.cpp.o"

# External object files for target generatecnatrees
generatecnatrees_EXTERNAL_OBJECTS =

generatecnatrees: CMakeFiles/generatecnatrees.dir/src/arg_parser.cpp.o
generatecnatrees: CMakeFiles/generatecnatrees.dir/src/generatecnatreesmain.cpp.o
generatecnatrees: CMakeFiles/generatecnatrees.dir/src/cnatree.cpp.o
generatecnatrees: CMakeFiles/generatecnatrees.dir/src/cnagraph.cpp.o
generatecnatrees: CMakeFiles/generatecnatrees.dir/src/genotypetree.cpp.o
generatecnatrees: CMakeFiles/generatecnatrees.dir/src/utils.cpp.o
generatecnatrees: CMakeFiles/generatecnatrees.dir/src/basematrix.cpp.o
generatecnatrees: CMakeFiles/generatecnatrees.dir/src/basetree.cpp.o
generatecnatrees: CMakeFiles/generatecnatrees.dir/build.make
generatecnatrees: /usr/lib64/libboost_thread-mt.so
generatecnatrees: /usr/lib64/libboost_system-mt.so
generatecnatrees: /usr/lib64/libboost_filesystem-mt.so
generatecnatrees: CMakeFiles/generatecnatrees.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable generatecnatrees"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/generatecnatrees.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/generatecnatrees.dir/build: generatecnatrees
.PHONY : CMakeFiles/generatecnatrees.dir/build

CMakeFiles/generatecnatrees.dir/requires: CMakeFiles/generatecnatrees.dir/src/arg_parser.cpp.o.requires
CMakeFiles/generatecnatrees.dir/requires: CMakeFiles/generatecnatrees.dir/src/generatecnatreesmain.cpp.o.requires
CMakeFiles/generatecnatrees.dir/requires: CMakeFiles/generatecnatrees.dir/src/cnatree.cpp.o.requires
CMakeFiles/generatecnatrees.dir/requires: CMakeFiles/generatecnatrees.dir/src/cnagraph.cpp.o.requires
CMakeFiles/generatecnatrees.dir/requires: CMakeFiles/generatecnatrees.dir/src/genotypetree.cpp.o.requires
CMakeFiles/generatecnatrees.dir/requires: CMakeFiles/generatecnatrees.dir/src/utils.cpp.o.requires
CMakeFiles/generatecnatrees.dir/requires: CMakeFiles/generatecnatrees.dir/src/basematrix.cpp.o.requires
CMakeFiles/generatecnatrees.dir/requires: CMakeFiles/generatecnatrees.dir/src/basetree.cpp.o.requires
.PHONY : CMakeFiles/generatecnatrees.dir/requires

CMakeFiles/generatecnatrees.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/generatecnatrees.dir/cmake_clean.cmake
.PHONY : CMakeFiles/generatecnatrees.dir/clean

CMakeFiles/generatecnatrees.dir/depend:
	cd /scratch/data/anna/clonesim/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /scratch/data/anna/clonesim /scratch/data/anna/clonesim /scratch/data/anna/clonesim/build /scratch/data/anna/clonesim/build /scratch/data/anna/clonesim/build/CMakeFiles/generatecnatrees.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/generatecnatrees.dir/depend

