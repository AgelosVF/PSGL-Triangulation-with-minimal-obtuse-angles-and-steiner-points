# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.31

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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/Aggelos/Projects/Project2/make

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/Aggelos/Projects/Project2/make

# Utility rule file for CGAL_Qt5_moc_and_resources_autogen.

# Include any custom commands dependencies for this target.
include CMakeFiles/CGAL_Qt5_moc_and_resources_autogen.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/CGAL_Qt5_moc_and_resources_autogen.dir/progress.make

CMakeFiles/CGAL_Qt5_moc_and_resources_autogen: CGAL_Qt5_moc_and_resources_autogen/timestamp

CGAL_Qt5_moc_and_resources_autogen/timestamp: /usr/bin/moc
CGAL_Qt5_moc_and_resources_autogen/timestamp: CMakeFiles/CGAL_Qt5_moc_and_resources_autogen.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --blue --bold --progress-dir=/home/Aggelos/Projects/Project2/make/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Automatic MOC for target CGAL_Qt5_moc_and_resources"
	/usr/bin/cmake -E cmake_autogen /home/Aggelos/Projects/Project2/make/CMakeFiles/CGAL_Qt5_moc_and_resources_autogen.dir/AutogenInfo.json ""
	/usr/bin/cmake -E touch /home/Aggelos/Projects/Project2/make/CGAL_Qt5_moc_and_resources_autogen/timestamp

CMakeFiles/CGAL_Qt5_moc_and_resources_autogen.dir/codegen:
.PHONY : CMakeFiles/CGAL_Qt5_moc_and_resources_autogen.dir/codegen

CGAL_Qt5_moc_and_resources_autogen: CGAL_Qt5_moc_and_resources_autogen/timestamp
CGAL_Qt5_moc_and_resources_autogen: CMakeFiles/CGAL_Qt5_moc_and_resources_autogen
CGAL_Qt5_moc_and_resources_autogen: CMakeFiles/CGAL_Qt5_moc_and_resources_autogen.dir/build.make
.PHONY : CGAL_Qt5_moc_and_resources_autogen

# Rule to build all files generated by this target.
CMakeFiles/CGAL_Qt5_moc_and_resources_autogen.dir/build: CGAL_Qt5_moc_and_resources_autogen
.PHONY : CMakeFiles/CGAL_Qt5_moc_and_resources_autogen.dir/build

CMakeFiles/CGAL_Qt5_moc_and_resources_autogen.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/CGAL_Qt5_moc_and_resources_autogen.dir/cmake_clean.cmake
.PHONY : CMakeFiles/CGAL_Qt5_moc_and_resources_autogen.dir/clean

CMakeFiles/CGAL_Qt5_moc_and_resources_autogen.dir/depend:
	cd /home/Aggelos/Projects/Project2/make && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/Aggelos/Projects/Project2/make /home/Aggelos/Projects/Project2/make /home/Aggelos/Projects/Project2/make /home/Aggelos/Projects/Project2/make /home/Aggelos/Projects/Project2/make/CMakeFiles/CGAL_Qt5_moc_and_resources_autogen.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : CMakeFiles/CGAL_Qt5_moc_and_resources_autogen.dir/depend

