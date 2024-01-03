# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.28

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
CMAKE_COMMAND = /usr/local/Homebrew/Cellar/cmake/3.28.0/bin/cmake

# The command to remove a file.
RM = /usr/local/Homebrew/Cellar/cmake/3.28.0/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/truonggiangdo/Library/CloudStorage/OneDrive-PolitechnikaWarszawska/Thesis/Survey-and-comparison-of-SSSP-algorithms/src/Dijkstra

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/truonggiangdo/Library/CloudStorage/OneDrive-PolitechnikaWarszawska/Thesis/Survey-and-comparison-of-SSSP-algorithms/src/Dijkstra/build

# Include any dependencies generated for this target.
include CMakeFiles/Dijkstra.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/Dijkstra.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/Dijkstra.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/Dijkstra.dir/flags.make

CMakeFiles/Dijkstra.dir/main.cpp.o: CMakeFiles/Dijkstra.dir/flags.make
CMakeFiles/Dijkstra.dir/main.cpp.o: /Users/truonggiangdo/Library/CloudStorage/OneDrive-PolitechnikaWarszawska/Thesis/Survey-and-comparison-of-SSSP-algorithms/src/Dijkstra/main.cpp
CMakeFiles/Dijkstra.dir/main.cpp.o: CMakeFiles/Dijkstra.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/truonggiangdo/Library/CloudStorage/OneDrive-PolitechnikaWarszawska/Thesis/Survey-and-comparison-of-SSSP-algorithms/src/Dijkstra/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/Dijkstra.dir/main.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/Dijkstra.dir/main.cpp.o -MF CMakeFiles/Dijkstra.dir/main.cpp.o.d -o CMakeFiles/Dijkstra.dir/main.cpp.o -c /Users/truonggiangdo/Library/CloudStorage/OneDrive-PolitechnikaWarszawska/Thesis/Survey-and-comparison-of-SSSP-algorithms/src/Dijkstra/main.cpp

CMakeFiles/Dijkstra.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/Dijkstra.dir/main.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/truonggiangdo/Library/CloudStorage/OneDrive-PolitechnikaWarszawska/Thesis/Survey-and-comparison-of-SSSP-algorithms/src/Dijkstra/main.cpp > CMakeFiles/Dijkstra.dir/main.cpp.i

CMakeFiles/Dijkstra.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/Dijkstra.dir/main.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/truonggiangdo/Library/CloudStorage/OneDrive-PolitechnikaWarszawska/Thesis/Survey-and-comparison-of-SSSP-algorithms/src/Dijkstra/main.cpp -o CMakeFiles/Dijkstra.dir/main.cpp.s

# Object files for target Dijkstra
Dijkstra_OBJECTS = \
"CMakeFiles/Dijkstra.dir/main.cpp.o"

# External object files for target Dijkstra
Dijkstra_EXTERNAL_OBJECTS =

Dijkstra: CMakeFiles/Dijkstra.dir/main.cpp.o
Dijkstra: CMakeFiles/Dijkstra.dir/build.make
Dijkstra: CMakeFiles/Dijkstra.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=/Users/truonggiangdo/Library/CloudStorage/OneDrive-PolitechnikaWarszawska/Thesis/Survey-and-comparison-of-SSSP-algorithms/src/Dijkstra/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable Dijkstra"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/Dijkstra.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/Dijkstra.dir/build: Dijkstra
.PHONY : CMakeFiles/Dijkstra.dir/build

CMakeFiles/Dijkstra.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/Dijkstra.dir/cmake_clean.cmake
.PHONY : CMakeFiles/Dijkstra.dir/clean

CMakeFiles/Dijkstra.dir/depend:
	cd /Users/truonggiangdo/Library/CloudStorage/OneDrive-PolitechnikaWarszawska/Thesis/Survey-and-comparison-of-SSSP-algorithms/src/Dijkstra/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/truonggiangdo/Library/CloudStorage/OneDrive-PolitechnikaWarszawska/Thesis/Survey-and-comparison-of-SSSP-algorithms/src/Dijkstra /Users/truonggiangdo/Library/CloudStorage/OneDrive-PolitechnikaWarszawska/Thesis/Survey-and-comparison-of-SSSP-algorithms/src/Dijkstra /Users/truonggiangdo/Library/CloudStorage/OneDrive-PolitechnikaWarszawska/Thesis/Survey-and-comparison-of-SSSP-algorithms/src/Dijkstra/build /Users/truonggiangdo/Library/CloudStorage/OneDrive-PolitechnikaWarszawska/Thesis/Survey-and-comparison-of-SSSP-algorithms/src/Dijkstra/build /Users/truonggiangdo/Library/CloudStorage/OneDrive-PolitechnikaWarszawska/Thesis/Survey-and-comparison-of-SSSP-algorithms/src/Dijkstra/build/CMakeFiles/Dijkstra.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : CMakeFiles/Dijkstra.dir/depend

