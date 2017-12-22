# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.7

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


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
CMAKE_COMMAND = /opt/local/bin/cmake

# The command to remove a file.
RM = /opt/local/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/angus1/Programs/MHDmpi/src

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/angus1/Programs/MHDmpi/physicsMods/mhdShock

# Include any dependencies generated for this target.
include CMakeFiles/mhd.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/mhd.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/mhd.dir/flags.make

CMakeFiles/mhd.dir/main.cpp.o: CMakeFiles/mhd.dir/flags.make
CMakeFiles/mhd.dir/main.cpp.o: /Users/angus1/Programs/MHDmpi/src/main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/angus1/Programs/MHDmpi/physicsMods/mhdShock/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/mhd.dir/main.cpp.o"
	mpicxx   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/mhd.dir/main.cpp.o -c /Users/angus1/Programs/MHDmpi/src/main.cpp

CMakeFiles/mhd.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mhd.dir/main.cpp.i"
	mpicxx  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/angus1/Programs/MHDmpi/src/main.cpp > CMakeFiles/mhd.dir/main.cpp.i

CMakeFiles/mhd.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mhd.dir/main.cpp.s"
	mpicxx  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/angus1/Programs/MHDmpi/src/main.cpp -o CMakeFiles/mhd.dir/main.cpp.s

CMakeFiles/mhd.dir/main.cpp.o.requires:

.PHONY : CMakeFiles/mhd.dir/main.cpp.o.requires

CMakeFiles/mhd.dir/main.cpp.o.provides: CMakeFiles/mhd.dir/main.cpp.o.requires
	$(MAKE) -f CMakeFiles/mhd.dir/build.make CMakeFiles/mhd.dir/main.cpp.o.provides.build
.PHONY : CMakeFiles/mhd.dir/main.cpp.o.provides

CMakeFiles/mhd.dir/main.cpp.o.provides.build: CMakeFiles/mhd.dir/main.cpp.o


CMakeFiles/mhd.dir/domainGrid.cpp.o: CMakeFiles/mhd.dir/flags.make
CMakeFiles/mhd.dir/domainGrid.cpp.o: /Users/angus1/Programs/MHDmpi/src/domainGrid.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/angus1/Programs/MHDmpi/physicsMods/mhdShock/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/mhd.dir/domainGrid.cpp.o"
	mpicxx   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/mhd.dir/domainGrid.cpp.o -c /Users/angus1/Programs/MHDmpi/src/domainGrid.cpp

CMakeFiles/mhd.dir/domainGrid.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mhd.dir/domainGrid.cpp.i"
	mpicxx  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/angus1/Programs/MHDmpi/src/domainGrid.cpp > CMakeFiles/mhd.dir/domainGrid.cpp.i

CMakeFiles/mhd.dir/domainGrid.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mhd.dir/domainGrid.cpp.s"
	mpicxx  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/angus1/Programs/MHDmpi/src/domainGrid.cpp -o CMakeFiles/mhd.dir/domainGrid.cpp.s

CMakeFiles/mhd.dir/domainGrid.cpp.o.requires:

.PHONY : CMakeFiles/mhd.dir/domainGrid.cpp.o.requires

CMakeFiles/mhd.dir/domainGrid.cpp.o.provides: CMakeFiles/mhd.dir/domainGrid.cpp.o.requires
	$(MAKE) -f CMakeFiles/mhd.dir/build.make CMakeFiles/mhd.dir/domainGrid.cpp.o.provides.build
.PHONY : CMakeFiles/mhd.dir/domainGrid.cpp.o.provides

CMakeFiles/mhd.dir/domainGrid.cpp.o.provides.build: CMakeFiles/mhd.dir/domainGrid.cpp.o


CMakeFiles/mhd.dir/vectorMath.cpp.o: CMakeFiles/mhd.dir/flags.make
CMakeFiles/mhd.dir/vectorMath.cpp.o: /Users/angus1/Programs/MHDmpi/src/vectorMath.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/angus1/Programs/MHDmpi/physicsMods/mhdShock/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/mhd.dir/vectorMath.cpp.o"
	mpicxx   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/mhd.dir/vectorMath.cpp.o -c /Users/angus1/Programs/MHDmpi/src/vectorMath.cpp

CMakeFiles/mhd.dir/vectorMath.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mhd.dir/vectorMath.cpp.i"
	mpicxx  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/angus1/Programs/MHDmpi/src/vectorMath.cpp > CMakeFiles/mhd.dir/vectorMath.cpp.i

CMakeFiles/mhd.dir/vectorMath.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mhd.dir/vectorMath.cpp.s"
	mpicxx  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/angus1/Programs/MHDmpi/src/vectorMath.cpp -o CMakeFiles/mhd.dir/vectorMath.cpp.s

CMakeFiles/mhd.dir/vectorMath.cpp.o.requires:

.PHONY : CMakeFiles/mhd.dir/vectorMath.cpp.o.requires

CMakeFiles/mhd.dir/vectorMath.cpp.o.provides: CMakeFiles/mhd.dir/vectorMath.cpp.o.requires
	$(MAKE) -f CMakeFiles/mhd.dir/build.make CMakeFiles/mhd.dir/vectorMath.cpp.o.provides.build
.PHONY : CMakeFiles/mhd.dir/vectorMath.cpp.o.provides

CMakeFiles/mhd.dir/vectorMath.cpp.o.provides.build: CMakeFiles/mhd.dir/vectorMath.cpp.o


CMakeFiles/mhd.dir/timeDomain.cpp.o: CMakeFiles/mhd.dir/flags.make
CMakeFiles/mhd.dir/timeDomain.cpp.o: /Users/angus1/Programs/MHDmpi/src/timeDomain.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/angus1/Programs/MHDmpi/physicsMods/mhdShock/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/mhd.dir/timeDomain.cpp.o"
	mpicxx   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/mhd.dir/timeDomain.cpp.o -c /Users/angus1/Programs/MHDmpi/src/timeDomain.cpp

CMakeFiles/mhd.dir/timeDomain.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mhd.dir/timeDomain.cpp.i"
	mpicxx  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/angus1/Programs/MHDmpi/src/timeDomain.cpp > CMakeFiles/mhd.dir/timeDomain.cpp.i

CMakeFiles/mhd.dir/timeDomain.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mhd.dir/timeDomain.cpp.s"
	mpicxx  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/angus1/Programs/MHDmpi/src/timeDomain.cpp -o CMakeFiles/mhd.dir/timeDomain.cpp.s

CMakeFiles/mhd.dir/timeDomain.cpp.o.requires:

.PHONY : CMakeFiles/mhd.dir/timeDomain.cpp.o.requires

CMakeFiles/mhd.dir/timeDomain.cpp.o.provides: CMakeFiles/mhd.dir/timeDomain.cpp.o.requires
	$(MAKE) -f CMakeFiles/mhd.dir/build.make CMakeFiles/mhd.dir/timeDomain.cpp.o.provides.build
.PHONY : CMakeFiles/mhd.dir/timeDomain.cpp.o.provides

CMakeFiles/mhd.dir/timeDomain.cpp.o.provides.build: CMakeFiles/mhd.dir/timeDomain.cpp.o


CMakeFiles/mhd.dir/HDF5dataFile.cpp.o: CMakeFiles/mhd.dir/flags.make
CMakeFiles/mhd.dir/HDF5dataFile.cpp.o: /Users/angus1/Programs/MHDmpi/src/HDF5dataFile.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/angus1/Programs/MHDmpi/physicsMods/mhdShock/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMakeFiles/mhd.dir/HDF5dataFile.cpp.o"
	mpicxx   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/mhd.dir/HDF5dataFile.cpp.o -c /Users/angus1/Programs/MHDmpi/src/HDF5dataFile.cpp

CMakeFiles/mhd.dir/HDF5dataFile.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mhd.dir/HDF5dataFile.cpp.i"
	mpicxx  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/angus1/Programs/MHDmpi/src/HDF5dataFile.cpp > CMakeFiles/mhd.dir/HDF5dataFile.cpp.i

CMakeFiles/mhd.dir/HDF5dataFile.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mhd.dir/HDF5dataFile.cpp.s"
	mpicxx  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/angus1/Programs/MHDmpi/src/HDF5dataFile.cpp -o CMakeFiles/mhd.dir/HDF5dataFile.cpp.s

CMakeFiles/mhd.dir/HDF5dataFile.cpp.o.requires:

.PHONY : CMakeFiles/mhd.dir/HDF5dataFile.cpp.o.requires

CMakeFiles/mhd.dir/HDF5dataFile.cpp.o.provides: CMakeFiles/mhd.dir/HDF5dataFile.cpp.o.requires
	$(MAKE) -f CMakeFiles/mhd.dir/build.make CMakeFiles/mhd.dir/HDF5dataFile.cpp.o.provides.build
.PHONY : CMakeFiles/mhd.dir/HDF5dataFile.cpp.o.provides

CMakeFiles/mhd.dir/HDF5dataFile.cpp.o.provides.build: CMakeFiles/mhd.dir/HDF5dataFile.cpp.o


CMakeFiles/mhd.dir/jsoncpp.cpp.o: CMakeFiles/mhd.dir/flags.make
CMakeFiles/mhd.dir/jsoncpp.cpp.o: /Users/angus1/Programs/MHDmpi/src/jsoncpp.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/angus1/Programs/MHDmpi/physicsMods/mhdShock/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object CMakeFiles/mhd.dir/jsoncpp.cpp.o"
	mpicxx   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/mhd.dir/jsoncpp.cpp.o -c /Users/angus1/Programs/MHDmpi/src/jsoncpp.cpp

CMakeFiles/mhd.dir/jsoncpp.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mhd.dir/jsoncpp.cpp.i"
	mpicxx  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/angus1/Programs/MHDmpi/src/jsoncpp.cpp > CMakeFiles/mhd.dir/jsoncpp.cpp.i

CMakeFiles/mhd.dir/jsoncpp.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mhd.dir/jsoncpp.cpp.s"
	mpicxx  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/angus1/Programs/MHDmpi/src/jsoncpp.cpp -o CMakeFiles/mhd.dir/jsoncpp.cpp.s

CMakeFiles/mhd.dir/jsoncpp.cpp.o.requires:

.PHONY : CMakeFiles/mhd.dir/jsoncpp.cpp.o.requires

CMakeFiles/mhd.dir/jsoncpp.cpp.o.provides: CMakeFiles/mhd.dir/jsoncpp.cpp.o.requires
	$(MAKE) -f CMakeFiles/mhd.dir/build.make CMakeFiles/mhd.dir/jsoncpp.cpp.o.provides.build
.PHONY : CMakeFiles/mhd.dir/jsoncpp.cpp.o.provides

CMakeFiles/mhd.dir/jsoncpp.cpp.o.provides.build: CMakeFiles/mhd.dir/jsoncpp.cpp.o


CMakeFiles/mhd.dir/mhdShock.cpp.o: CMakeFiles/mhd.dir/flags.make
CMakeFiles/mhd.dir/mhdShock.cpp.o: mhdShock.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/angus1/Programs/MHDmpi/physicsMods/mhdShock/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building CXX object CMakeFiles/mhd.dir/mhdShock.cpp.o"
	mpicxx   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/mhd.dir/mhdShock.cpp.o -c /Users/angus1/Programs/MHDmpi/physicsMods/mhdShock/mhdShock.cpp

CMakeFiles/mhd.dir/mhdShock.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mhd.dir/mhdShock.cpp.i"
	mpicxx  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/angus1/Programs/MHDmpi/physicsMods/mhdShock/mhdShock.cpp > CMakeFiles/mhd.dir/mhdShock.cpp.i

CMakeFiles/mhd.dir/mhdShock.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mhd.dir/mhdShock.cpp.s"
	mpicxx  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/angus1/Programs/MHDmpi/physicsMods/mhdShock/mhdShock.cpp -o CMakeFiles/mhd.dir/mhdShock.cpp.s

CMakeFiles/mhd.dir/mhdShock.cpp.o.requires:

.PHONY : CMakeFiles/mhd.dir/mhdShock.cpp.o.requires

CMakeFiles/mhd.dir/mhdShock.cpp.o.provides: CMakeFiles/mhd.dir/mhdShock.cpp.o.requires
	$(MAKE) -f CMakeFiles/mhd.dir/build.make CMakeFiles/mhd.dir/mhdShock.cpp.o.provides.build
.PHONY : CMakeFiles/mhd.dir/mhdShock.cpp.o.provides

CMakeFiles/mhd.dir/mhdShock.cpp.o.provides.build: CMakeFiles/mhd.dir/mhdShock.cpp.o


# Object files for target mhd
mhd_OBJECTS = \
"CMakeFiles/mhd.dir/main.cpp.o" \
"CMakeFiles/mhd.dir/domainGrid.cpp.o" \
"CMakeFiles/mhd.dir/vectorMath.cpp.o" \
"CMakeFiles/mhd.dir/timeDomain.cpp.o" \
"CMakeFiles/mhd.dir/HDF5dataFile.cpp.o" \
"CMakeFiles/mhd.dir/jsoncpp.cpp.o" \
"CMakeFiles/mhd.dir/mhdShock.cpp.o"

# External object files for target mhd
mhd_EXTERNAL_OBJECTS =

mhd: CMakeFiles/mhd.dir/main.cpp.o
mhd: CMakeFiles/mhd.dir/domainGrid.cpp.o
mhd: CMakeFiles/mhd.dir/vectorMath.cpp.o
mhd: CMakeFiles/mhd.dir/timeDomain.cpp.o
mhd: CMakeFiles/mhd.dir/HDF5dataFile.cpp.o
mhd: CMakeFiles/mhd.dir/jsoncpp.cpp.o
mhd: CMakeFiles/mhd.dir/mhdShock.cpp.o
mhd: CMakeFiles/mhd.dir/build.make
mhd: /opt/local/lib/openmpi-gcc5/libmpi_cxx.dylib
mhd: /opt/local/lib/openmpi-gcc5/libmpi.dylib
mhd: CMakeFiles/mhd.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/angus1/Programs/MHDmpi/physicsMods/mhdShock/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Linking CXX executable mhd"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/mhd.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/mhd.dir/build: mhd

.PHONY : CMakeFiles/mhd.dir/build

CMakeFiles/mhd.dir/requires: CMakeFiles/mhd.dir/main.cpp.o.requires
CMakeFiles/mhd.dir/requires: CMakeFiles/mhd.dir/domainGrid.cpp.o.requires
CMakeFiles/mhd.dir/requires: CMakeFiles/mhd.dir/vectorMath.cpp.o.requires
CMakeFiles/mhd.dir/requires: CMakeFiles/mhd.dir/timeDomain.cpp.o.requires
CMakeFiles/mhd.dir/requires: CMakeFiles/mhd.dir/HDF5dataFile.cpp.o.requires
CMakeFiles/mhd.dir/requires: CMakeFiles/mhd.dir/jsoncpp.cpp.o.requires
CMakeFiles/mhd.dir/requires: CMakeFiles/mhd.dir/mhdShock.cpp.o.requires

.PHONY : CMakeFiles/mhd.dir/requires

CMakeFiles/mhd.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/mhd.dir/cmake_clean.cmake
.PHONY : CMakeFiles/mhd.dir/clean

CMakeFiles/mhd.dir/depend:
	cd /Users/angus1/Programs/MHDmpi/physicsMods/mhdShock && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/angus1/Programs/MHDmpi/src /Users/angus1/Programs/MHDmpi/src /Users/angus1/Programs/MHDmpi/physicsMods/mhdShock /Users/angus1/Programs/MHDmpi/physicsMods/mhdShock /Users/angus1/Programs/MHDmpi/physicsMods/mhdShock/CMakeFiles/mhd.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/mhd.dir/depend

