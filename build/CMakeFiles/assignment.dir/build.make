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
CMAKE_SOURCE_DIR = /v/filer4b/v38q001/rescobar/graphics/assignment4

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /v/filer4b/v38q001/rescobar/graphics/assignment4/build

# Include any dependencies generated for this target.
include CMakeFiles/assignment.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/assignment.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/assignment.dir/flags.make

CMakeFiles/assignment.dir/assignment.cc.o: CMakeFiles/assignment.dir/flags.make
CMakeFiles/assignment.dir/assignment.cc.o: ../assignment.cc
	$(CMAKE_COMMAND) -E cmake_progress_report /v/filer4b/v38q001/rescobar/graphics/assignment4/build/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/assignment.dir/assignment.cc.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/assignment.dir/assignment.cc.o -c /v/filer4b/v38q001/rescobar/graphics/assignment4/assignment.cc

CMakeFiles/assignment.dir/assignment.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/assignment.dir/assignment.cc.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /v/filer4b/v38q001/rescobar/graphics/assignment4/assignment.cc > CMakeFiles/assignment.dir/assignment.cc.i

CMakeFiles/assignment.dir/assignment.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/assignment.dir/assignment.cc.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /v/filer4b/v38q001/rescobar/graphics/assignment4/assignment.cc -o CMakeFiles/assignment.dir/assignment.cc.s

CMakeFiles/assignment.dir/assignment.cc.o.requires:
.PHONY : CMakeFiles/assignment.dir/assignment.cc.o.requires

CMakeFiles/assignment.dir/assignment.cc.o.provides: CMakeFiles/assignment.dir/assignment.cc.o.requires
	$(MAKE) -f CMakeFiles/assignment.dir/build.make CMakeFiles/assignment.dir/assignment.cc.o.provides.build
.PHONY : CMakeFiles/assignment.dir/assignment.cc.o.provides

CMakeFiles/assignment.dir/assignment.cc.o.provides.build: CMakeFiles/assignment.dir/assignment.cc.o

# Object files for target assignment
assignment_OBJECTS = \
"CMakeFiles/assignment.dir/assignment.cc.o"

# External object files for target assignment
assignment_EXTERNAL_OBJECTS =

bin/assignment: CMakeFiles/assignment.dir/assignment.cc.o
bin/assignment: CMakeFiles/assignment.dir/build.make
bin/assignment: /usr/lib/x86_64-linux-gnu/libGLEW.so
bin/assignment: /usr/lib/x86_64-linux-gnu/libGL.so
bin/assignment: /usr/lib/x86_64-linux-gnu/libGLEW.so
bin/assignment: /usr/lib/x86_64-linux-gnu/libjpeg.so
bin/assignment: /usr/lib/x86_64-linux-gnu/libGL.so
bin/assignment: /usr/lib/x86_64-linux-gnu/libjpeg.so
bin/assignment: CMakeFiles/assignment.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable bin/assignment"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/assignment.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/assignment.dir/build: bin/assignment
.PHONY : CMakeFiles/assignment.dir/build

CMakeFiles/assignment.dir/requires: CMakeFiles/assignment.dir/assignment.cc.o.requires
.PHONY : CMakeFiles/assignment.dir/requires

CMakeFiles/assignment.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/assignment.dir/cmake_clean.cmake
.PHONY : CMakeFiles/assignment.dir/clean

CMakeFiles/assignment.dir/depend:
	cd /v/filer4b/v38q001/rescobar/graphics/assignment4/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /v/filer4b/v38q001/rescobar/graphics/assignment4 /v/filer4b/v38q001/rescobar/graphics/assignment4 /v/filer4b/v38q001/rescobar/graphics/assignment4/build /v/filer4b/v38q001/rescobar/graphics/assignment4/build /v/filer4b/v38q001/rescobar/graphics/assignment4/build/CMakeFiles/assignment.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/assignment.dir/depend

