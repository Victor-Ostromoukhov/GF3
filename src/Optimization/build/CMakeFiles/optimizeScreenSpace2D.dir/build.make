# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.17

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
CMAKE_COMMAND = /usr/local/Cellar/cmake/3.17.2/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/3.17.2/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/ostrom/SobolPlusPlus/src/Optimization

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/ostrom/SobolPlusPlus/src/Optimization/build

# Include any dependencies generated for this target.
include CMakeFiles/optimizeScreenSpace2D.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/optimizeScreenSpace2D.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/optimizeScreenSpace2D.dir/flags.make

CMakeFiles/optimizeScreenSpace2D.dir/optimizeScreenSpace2D.cpp.o: CMakeFiles/optimizeScreenSpace2D.dir/flags.make
CMakeFiles/optimizeScreenSpace2D.dir/optimizeScreenSpace2D.cpp.o: ../optimizeScreenSpace2D.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/ostrom/SobolPlusPlus/src/Optimization/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/optimizeScreenSpace2D.dir/optimizeScreenSpace2D.cpp.o"
	/Applications/Xcode.app/Contents/Developer/usr/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/optimizeScreenSpace2D.dir/optimizeScreenSpace2D.cpp.o -c /Users/ostrom/SobolPlusPlus/src/Optimization/optimizeScreenSpace2D.cpp

CMakeFiles/optimizeScreenSpace2D.dir/optimizeScreenSpace2D.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/optimizeScreenSpace2D.dir/optimizeScreenSpace2D.cpp.i"
	/Applications/Xcode.app/Contents/Developer/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/ostrom/SobolPlusPlus/src/Optimization/optimizeScreenSpace2D.cpp > CMakeFiles/optimizeScreenSpace2D.dir/optimizeScreenSpace2D.cpp.i

CMakeFiles/optimizeScreenSpace2D.dir/optimizeScreenSpace2D.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/optimizeScreenSpace2D.dir/optimizeScreenSpace2D.cpp.s"
	/Applications/Xcode.app/Contents/Developer/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/ostrom/SobolPlusPlus/src/Optimization/optimizeScreenSpace2D.cpp -o CMakeFiles/optimizeScreenSpace2D.dir/optimizeScreenSpace2D.cpp.s

# Object files for target optimizeScreenSpace2D
optimizeScreenSpace2D_OBJECTS = \
"CMakeFiles/optimizeScreenSpace2D.dir/optimizeScreenSpace2D.cpp.o"

# External object files for target optimizeScreenSpace2D
optimizeScreenSpace2D_EXTERNAL_OBJECTS =

optimizeScreenSpace2D: CMakeFiles/optimizeScreenSpace2D.dir/optimizeScreenSpace2D.cpp.o
optimizeScreenSpace2D: CMakeFiles/optimizeScreenSpace2D.dir/build.make
optimizeScreenSpace2D: CMakeFiles/optimizeScreenSpace2D.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/ostrom/SobolPlusPlus/src/Optimization/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable optimizeScreenSpace2D"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/optimizeScreenSpace2D.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/optimizeScreenSpace2D.dir/build: optimizeScreenSpace2D

.PHONY : CMakeFiles/optimizeScreenSpace2D.dir/build

CMakeFiles/optimizeScreenSpace2D.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/optimizeScreenSpace2D.dir/cmake_clean.cmake
.PHONY : CMakeFiles/optimizeScreenSpace2D.dir/clean

CMakeFiles/optimizeScreenSpace2D.dir/depend:
	cd /Users/ostrom/SobolPlusPlus/src/Optimization/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/ostrom/SobolPlusPlus/src/Optimization /Users/ostrom/SobolPlusPlus/src/Optimization /Users/ostrom/SobolPlusPlus/src/Optimization/build /Users/ostrom/SobolPlusPlus/src/Optimization/build /Users/ostrom/SobolPlusPlus/src/Optimization/build/CMakeFiles/optimizeScreenSpace2D.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/optimizeScreenSpace2D.dir/depend

