# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.19

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
CMAKE_COMMAND = /usr/local/Cellar/cmake/3.19.6/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/3.19.6/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/Daniel/Documents/CODE/Research/Kansas/Cycler/C++

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/Daniel/Documents/CODE/Research/Kansas/Cycler/C++/build

# Include any dependencies generated for this target.
include CMakeFiles/TestDebug.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/TestDebug.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/TestDebug.dir/flags.make

CMakeFiles/TestDebug.dir/src/LambertTest.cpp.o: CMakeFiles/TestDebug.dir/flags.make
CMakeFiles/TestDebug.dir/src/LambertTest.cpp.o: ../src/LambertTest.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/Daniel/Documents/CODE/Research/Kansas/Cycler/C++/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/TestDebug.dir/src/LambertTest.cpp.o"
	/usr/local/opt/llvm/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/TestDebug.dir/src/LambertTest.cpp.o -c /Users/Daniel/Documents/CODE/Research/Kansas/Cycler/C++/src/LambertTest.cpp

CMakeFiles/TestDebug.dir/src/LambertTest.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/TestDebug.dir/src/LambertTest.cpp.i"
	/usr/local/opt/llvm/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/Daniel/Documents/CODE/Research/Kansas/Cycler/C++/src/LambertTest.cpp > CMakeFiles/TestDebug.dir/src/LambertTest.cpp.i

CMakeFiles/TestDebug.dir/src/LambertTest.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/TestDebug.dir/src/LambertTest.cpp.s"
	/usr/local/opt/llvm/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/Daniel/Documents/CODE/Research/Kansas/Cycler/C++/src/LambertTest.cpp -o CMakeFiles/TestDebug.dir/src/LambertTest.cpp.s

CMakeFiles/TestDebug.dir/src/Lambert_Multi.cpp.o: CMakeFiles/TestDebug.dir/flags.make
CMakeFiles/TestDebug.dir/src/Lambert_Multi.cpp.o: ../src/Lambert_Multi.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/Daniel/Documents/CODE/Research/Kansas/Cycler/C++/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/TestDebug.dir/src/Lambert_Multi.cpp.o"
	/usr/local/opt/llvm/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/TestDebug.dir/src/Lambert_Multi.cpp.o -c /Users/Daniel/Documents/CODE/Research/Kansas/Cycler/C++/src/Lambert_Multi.cpp

CMakeFiles/TestDebug.dir/src/Lambert_Multi.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/TestDebug.dir/src/Lambert_Multi.cpp.i"
	/usr/local/opt/llvm/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/Daniel/Documents/CODE/Research/Kansas/Cycler/C++/src/Lambert_Multi.cpp > CMakeFiles/TestDebug.dir/src/Lambert_Multi.cpp.i

CMakeFiles/TestDebug.dir/src/Lambert_Multi.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/TestDebug.dir/src/Lambert_Multi.cpp.s"
	/usr/local/opt/llvm/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/Daniel/Documents/CODE/Research/Kansas/Cycler/C++/src/Lambert_Multi.cpp -o CMakeFiles/TestDebug.dir/src/Lambert_Multi.cpp.s

CMakeFiles/TestDebug.dir/src/Lambert_Battin.cpp.o: CMakeFiles/TestDebug.dir/flags.make
CMakeFiles/TestDebug.dir/src/Lambert_Battin.cpp.o: ../src/Lambert_Battin.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/Daniel/Documents/CODE/Research/Kansas/Cycler/C++/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/TestDebug.dir/src/Lambert_Battin.cpp.o"
	/usr/local/opt/llvm/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/TestDebug.dir/src/Lambert_Battin.cpp.o -c /Users/Daniel/Documents/CODE/Research/Kansas/Cycler/C++/src/Lambert_Battin.cpp

CMakeFiles/TestDebug.dir/src/Lambert_Battin.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/TestDebug.dir/src/Lambert_Battin.cpp.i"
	/usr/local/opt/llvm/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/Daniel/Documents/CODE/Research/Kansas/Cycler/C++/src/Lambert_Battin.cpp > CMakeFiles/TestDebug.dir/src/Lambert_Battin.cpp.i

CMakeFiles/TestDebug.dir/src/Lambert_Battin.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/TestDebug.dir/src/Lambert_Battin.cpp.s"
	/usr/local/opt/llvm/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/Daniel/Documents/CODE/Research/Kansas/Cycler/C++/src/Lambert_Battin.cpp -o CMakeFiles/TestDebug.dir/src/Lambert_Battin.cpp.s

CMakeFiles/TestDebug.dir/src/Lambert_Battin_Multi.cpp.o: CMakeFiles/TestDebug.dir/flags.make
CMakeFiles/TestDebug.dir/src/Lambert_Battin_Multi.cpp.o: ../src/Lambert_Battin_Multi.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/Daniel/Documents/CODE/Research/Kansas/Cycler/C++/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/TestDebug.dir/src/Lambert_Battin_Multi.cpp.o"
	/usr/local/opt/llvm/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/TestDebug.dir/src/Lambert_Battin_Multi.cpp.o -c /Users/Daniel/Documents/CODE/Research/Kansas/Cycler/C++/src/Lambert_Battin_Multi.cpp

CMakeFiles/TestDebug.dir/src/Lambert_Battin_Multi.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/TestDebug.dir/src/Lambert_Battin_Multi.cpp.i"
	/usr/local/opt/llvm/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/Daniel/Documents/CODE/Research/Kansas/Cycler/C++/src/Lambert_Battin_Multi.cpp > CMakeFiles/TestDebug.dir/src/Lambert_Battin_Multi.cpp.i

CMakeFiles/TestDebug.dir/src/Lambert_Battin_Multi.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/TestDebug.dir/src/Lambert_Battin_Multi.cpp.s"
	/usr/local/opt/llvm/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/Daniel/Documents/CODE/Research/Kansas/Cycler/C++/src/Lambert_Battin_Multi.cpp -o CMakeFiles/TestDebug.dir/src/Lambert_Battin_Multi.cpp.s

CMakeFiles/TestDebug.dir/src/helperFuncs.cpp.o: CMakeFiles/TestDebug.dir/flags.make
CMakeFiles/TestDebug.dir/src/helperFuncs.cpp.o: ../src/helperFuncs.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/Daniel/Documents/CODE/Research/Kansas/Cycler/C++/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMakeFiles/TestDebug.dir/src/helperFuncs.cpp.o"
	/usr/local/opt/llvm/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/TestDebug.dir/src/helperFuncs.cpp.o -c /Users/Daniel/Documents/CODE/Research/Kansas/Cycler/C++/src/helperFuncs.cpp

CMakeFiles/TestDebug.dir/src/helperFuncs.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/TestDebug.dir/src/helperFuncs.cpp.i"
	/usr/local/opt/llvm/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/Daniel/Documents/CODE/Research/Kansas/Cycler/C++/src/helperFuncs.cpp > CMakeFiles/TestDebug.dir/src/helperFuncs.cpp.i

CMakeFiles/TestDebug.dir/src/helperFuncs.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/TestDebug.dir/src/helperFuncs.cpp.s"
	/usr/local/opt/llvm/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/Daniel/Documents/CODE/Research/Kansas/Cycler/C++/src/helperFuncs.cpp -o CMakeFiles/TestDebug.dir/src/helperFuncs.cpp.s

# Object files for target TestDebug
TestDebug_OBJECTS = \
"CMakeFiles/TestDebug.dir/src/LambertTest.cpp.o" \
"CMakeFiles/TestDebug.dir/src/Lambert_Multi.cpp.o" \
"CMakeFiles/TestDebug.dir/src/Lambert_Battin.cpp.o" \
"CMakeFiles/TestDebug.dir/src/Lambert_Battin_Multi.cpp.o" \
"CMakeFiles/TestDebug.dir/src/helperFuncs.cpp.o"

# External object files for target TestDebug
TestDebug_EXTERNAL_OBJECTS =

TestDebug: CMakeFiles/TestDebug.dir/src/LambertTest.cpp.o
TestDebug: CMakeFiles/TestDebug.dir/src/Lambert_Multi.cpp.o
TestDebug: CMakeFiles/TestDebug.dir/src/Lambert_Battin.cpp.o
TestDebug: CMakeFiles/TestDebug.dir/src/Lambert_Battin_Multi.cpp.o
TestDebug: CMakeFiles/TestDebug.dir/src/helperFuncs.cpp.o
TestDebug: CMakeFiles/TestDebug.dir/build.make
TestDebug: CMakeFiles/TestDebug.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/Daniel/Documents/CODE/Research/Kansas/Cycler/C++/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Linking CXX executable TestDebug"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/TestDebug.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/TestDebug.dir/build: TestDebug

.PHONY : CMakeFiles/TestDebug.dir/build

CMakeFiles/TestDebug.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/TestDebug.dir/cmake_clean.cmake
.PHONY : CMakeFiles/TestDebug.dir/clean

CMakeFiles/TestDebug.dir/depend:
	cd /Users/Daniel/Documents/CODE/Research/Kansas/Cycler/C++/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/Daniel/Documents/CODE/Research/Kansas/Cycler/C++ /Users/Daniel/Documents/CODE/Research/Kansas/Cycler/C++ /Users/Daniel/Documents/CODE/Research/Kansas/Cycler/C++/build /Users/Daniel/Documents/CODE/Research/Kansas/Cycler/C++/build /Users/Daniel/Documents/CODE/Research/Kansas/Cycler/C++/build/CMakeFiles/TestDebug.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/TestDebug.dir/depend
