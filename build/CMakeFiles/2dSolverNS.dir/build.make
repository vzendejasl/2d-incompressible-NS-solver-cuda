# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/vzendejasl/Documents/CaltechCourses/CS179/CS179Project/232aCFDCode_CUDA

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/vzendejasl/Documents/CaltechCourses/CS179/CS179Project/232aCFDCode_CUDA/build

# Include any dependencies generated for this target.
include CMakeFiles/2dSolverNS.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/2dSolverNS.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/2dSolverNS.dir/flags.make

CMakeFiles/2dSolverNS.dir/main.cpp.o: CMakeFiles/2dSolverNS.dir/flags.make
CMakeFiles/2dSolverNS.dir/main.cpp.o: ../main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/vzendejasl/Documents/CaltechCourses/CS179/CS179Project/232aCFDCode_CUDA/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/2dSolverNS.dir/main.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/2dSolverNS.dir/main.cpp.o -c /home/vzendejasl/Documents/CaltechCourses/CS179/CS179Project/232aCFDCode_CUDA/main.cpp

CMakeFiles/2dSolverNS.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/2dSolverNS.dir/main.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/vzendejasl/Documents/CaltechCourses/CS179/CS179Project/232aCFDCode_CUDA/main.cpp > CMakeFiles/2dSolverNS.dir/main.cpp.i

CMakeFiles/2dSolverNS.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/2dSolverNS.dir/main.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/vzendejasl/Documents/CaltechCourses/CS179/CS179Project/232aCFDCode_CUDA/main.cpp -o CMakeFiles/2dSolverNS.dir/main.cpp.s

CMakeFiles/2dSolverNS.dir/src/compute_RHS.cpp.o: CMakeFiles/2dSolverNS.dir/flags.make
CMakeFiles/2dSolverNS.dir/src/compute_RHS.cpp.o: ../src/compute_RHS.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/vzendejasl/Documents/CaltechCourses/CS179/CS179Project/232aCFDCode_CUDA/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/2dSolverNS.dir/src/compute_RHS.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/2dSolverNS.dir/src/compute_RHS.cpp.o -c /home/vzendejasl/Documents/CaltechCourses/CS179/CS179Project/232aCFDCode_CUDA/src/compute_RHS.cpp

CMakeFiles/2dSolverNS.dir/src/compute_RHS.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/2dSolverNS.dir/src/compute_RHS.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/vzendejasl/Documents/CaltechCourses/CS179/CS179Project/232aCFDCode_CUDA/src/compute_RHS.cpp > CMakeFiles/2dSolverNS.dir/src/compute_RHS.cpp.i

CMakeFiles/2dSolverNS.dir/src/compute_RHS.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/2dSolverNS.dir/src/compute_RHS.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/vzendejasl/Documents/CaltechCourses/CS179/CS179Project/232aCFDCode_CUDA/src/compute_RHS.cpp -o CMakeFiles/2dSolverNS.dir/src/compute_RHS.cpp.s

CMakeFiles/2dSolverNS.dir/src/compute_curl.cpp.o: CMakeFiles/2dSolverNS.dir/flags.make
CMakeFiles/2dSolverNS.dir/src/compute_curl.cpp.o: ../src/compute_curl.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/vzendejasl/Documents/CaltechCourses/CS179/CS179Project/232aCFDCode_CUDA/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/2dSolverNS.dir/src/compute_curl.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/2dSolverNS.dir/src/compute_curl.cpp.o -c /home/vzendejasl/Documents/CaltechCourses/CS179/CS179Project/232aCFDCode_CUDA/src/compute_curl.cpp

CMakeFiles/2dSolverNS.dir/src/compute_curl.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/2dSolverNS.dir/src/compute_curl.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/vzendejasl/Documents/CaltechCourses/CS179/CS179Project/232aCFDCode_CUDA/src/compute_curl.cpp > CMakeFiles/2dSolverNS.dir/src/compute_curl.cpp.i

CMakeFiles/2dSolverNS.dir/src/compute_curl.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/2dSolverNS.dir/src/compute_curl.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/vzendejasl/Documents/CaltechCourses/CS179/CS179Project/232aCFDCode_CUDA/src/compute_curl.cpp -o CMakeFiles/2dSolverNS.dir/src/compute_curl.cpp.s

CMakeFiles/2dSolverNS.dir/src/compute_projection_step.cpp.o: CMakeFiles/2dSolverNS.dir/flags.make
CMakeFiles/2dSolverNS.dir/src/compute_projection_step.cpp.o: ../src/compute_projection_step.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/vzendejasl/Documents/CaltechCourses/CS179/CS179Project/232aCFDCode_CUDA/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/2dSolverNS.dir/src/compute_projection_step.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/2dSolverNS.dir/src/compute_projection_step.cpp.o -c /home/vzendejasl/Documents/CaltechCourses/CS179/CS179Project/232aCFDCode_CUDA/src/compute_projection_step.cpp

CMakeFiles/2dSolverNS.dir/src/compute_projection_step.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/2dSolverNS.dir/src/compute_projection_step.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/vzendejasl/Documents/CaltechCourses/CS179/CS179Project/232aCFDCode_CUDA/src/compute_projection_step.cpp > CMakeFiles/2dSolverNS.dir/src/compute_projection_step.cpp.i

CMakeFiles/2dSolverNS.dir/src/compute_projection_step.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/2dSolverNS.dir/src/compute_projection_step.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/vzendejasl/Documents/CaltechCourses/CS179/CS179Project/232aCFDCode_CUDA/src/compute_projection_step.cpp -o CMakeFiles/2dSolverNS.dir/src/compute_projection_step.cpp.s

CMakeFiles/2dSolverNS.dir/src/compute_RHS.cu.o: CMakeFiles/2dSolverNS.dir/flags.make
CMakeFiles/2dSolverNS.dir/src/compute_RHS.cu.o: ../src/compute_RHS.cu
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/vzendejasl/Documents/CaltechCourses/CS179/CS179Project/232aCFDCode_CUDA/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CUDA object CMakeFiles/2dSolverNS.dir/src/compute_RHS.cu.o"
	/usr/bin/nvcc  $(CUDA_DEFINES) $(CUDA_INCLUDES) $(CUDA_FLAGS) -x cu -c /home/vzendejasl/Documents/CaltechCourses/CS179/CS179Project/232aCFDCode_CUDA/src/compute_RHS.cu -o CMakeFiles/2dSolverNS.dir/src/compute_RHS.cu.o

CMakeFiles/2dSolverNS.dir/src/compute_RHS.cu.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CUDA source to CMakeFiles/2dSolverNS.dir/src/compute_RHS.cu.i"
	$(CMAKE_COMMAND) -E cmake_unimplemented_variable CMAKE_CUDA_CREATE_PREPROCESSED_SOURCE

CMakeFiles/2dSolverNS.dir/src/compute_RHS.cu.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CUDA source to assembly CMakeFiles/2dSolverNS.dir/src/compute_RHS.cu.s"
	$(CMAKE_COMMAND) -E cmake_unimplemented_variable CMAKE_CUDA_CREATE_ASSEMBLY_SOURCE

CMakeFiles/2dSolverNS.dir/src/compute_projection_step.cu.o: CMakeFiles/2dSolverNS.dir/flags.make
CMakeFiles/2dSolverNS.dir/src/compute_projection_step.cu.o: ../src/compute_projection_step.cu
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/vzendejasl/Documents/CaltechCourses/CS179/CS179Project/232aCFDCode_CUDA/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CUDA object CMakeFiles/2dSolverNS.dir/src/compute_projection_step.cu.o"
	/usr/bin/nvcc  $(CUDA_DEFINES) $(CUDA_INCLUDES) $(CUDA_FLAGS) -x cu -c /home/vzendejasl/Documents/CaltechCourses/CS179/CS179Project/232aCFDCode_CUDA/src/compute_projection_step.cu -o CMakeFiles/2dSolverNS.dir/src/compute_projection_step.cu.o

CMakeFiles/2dSolverNS.dir/src/compute_projection_step.cu.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CUDA source to CMakeFiles/2dSolverNS.dir/src/compute_projection_step.cu.i"
	$(CMAKE_COMMAND) -E cmake_unimplemented_variable CMAKE_CUDA_CREATE_PREPROCESSED_SOURCE

CMakeFiles/2dSolverNS.dir/src/compute_projection_step.cu.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CUDA source to assembly CMakeFiles/2dSolverNS.dir/src/compute_projection_step.cu.s"
	$(CMAKE_COMMAND) -E cmake_unimplemented_variable CMAKE_CUDA_CREATE_ASSEMBLY_SOURCE

CMakeFiles/2dSolverNS.dir/src/time_advance_RK3.cu.o: CMakeFiles/2dSolverNS.dir/flags.make
CMakeFiles/2dSolverNS.dir/src/time_advance_RK3.cu.o: ../src/time_advance_RK3.cu
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/vzendejasl/Documents/CaltechCourses/CS179/CS179Project/232aCFDCode_CUDA/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building CUDA object CMakeFiles/2dSolverNS.dir/src/time_advance_RK3.cu.o"
	/usr/bin/nvcc  $(CUDA_DEFINES) $(CUDA_INCLUDES) $(CUDA_FLAGS) -x cu -c /home/vzendejasl/Documents/CaltechCourses/CS179/CS179Project/232aCFDCode_CUDA/src/time_advance_RK3.cu -o CMakeFiles/2dSolverNS.dir/src/time_advance_RK3.cu.o

CMakeFiles/2dSolverNS.dir/src/time_advance_RK3.cu.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CUDA source to CMakeFiles/2dSolverNS.dir/src/time_advance_RK3.cu.i"
	$(CMAKE_COMMAND) -E cmake_unimplemented_variable CMAKE_CUDA_CREATE_PREPROCESSED_SOURCE

CMakeFiles/2dSolverNS.dir/src/time_advance_RK3.cu.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CUDA source to assembly CMakeFiles/2dSolverNS.dir/src/time_advance_RK3.cu.s"
	$(CMAKE_COMMAND) -E cmake_unimplemented_variable CMAKE_CUDA_CREATE_ASSEMBLY_SOURCE

CMakeFiles/2dSolverNS.dir/src/helper_functions.cpp.o: CMakeFiles/2dSolverNS.dir/flags.make
CMakeFiles/2dSolverNS.dir/src/helper_functions.cpp.o: ../src/helper_functions.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/vzendejasl/Documents/CaltechCourses/CS179/CS179Project/232aCFDCode_CUDA/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Building CXX object CMakeFiles/2dSolverNS.dir/src/helper_functions.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/2dSolverNS.dir/src/helper_functions.cpp.o -c /home/vzendejasl/Documents/CaltechCourses/CS179/CS179Project/232aCFDCode_CUDA/src/helper_functions.cpp

CMakeFiles/2dSolverNS.dir/src/helper_functions.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/2dSolverNS.dir/src/helper_functions.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/vzendejasl/Documents/CaltechCourses/CS179/CS179Project/232aCFDCode_CUDA/src/helper_functions.cpp > CMakeFiles/2dSolverNS.dir/src/helper_functions.cpp.i

CMakeFiles/2dSolverNS.dir/src/helper_functions.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/2dSolverNS.dir/src/helper_functions.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/vzendejasl/Documents/CaltechCourses/CS179/CS179Project/232aCFDCode_CUDA/src/helper_functions.cpp -o CMakeFiles/2dSolverNS.dir/src/helper_functions.cpp.s

CMakeFiles/2dSolverNS.dir/src/time_advance_RK3.cpp.o: CMakeFiles/2dSolverNS.dir/flags.make
CMakeFiles/2dSolverNS.dir/src/time_advance_RK3.cpp.o: ../src/time_advance_RK3.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/vzendejasl/Documents/CaltechCourses/CS179/CS179Project/232aCFDCode_CUDA/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_9) "Building CXX object CMakeFiles/2dSolverNS.dir/src/time_advance_RK3.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/2dSolverNS.dir/src/time_advance_RK3.cpp.o -c /home/vzendejasl/Documents/CaltechCourses/CS179/CS179Project/232aCFDCode_CUDA/src/time_advance_RK3.cpp

CMakeFiles/2dSolverNS.dir/src/time_advance_RK3.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/2dSolverNS.dir/src/time_advance_RK3.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/vzendejasl/Documents/CaltechCourses/CS179/CS179Project/232aCFDCode_CUDA/src/time_advance_RK3.cpp > CMakeFiles/2dSolverNS.dir/src/time_advance_RK3.cpp.i

CMakeFiles/2dSolverNS.dir/src/time_advance_RK3.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/2dSolverNS.dir/src/time_advance_RK3.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/vzendejasl/Documents/CaltechCourses/CS179/CS179Project/232aCFDCode_CUDA/src/time_advance_RK3.cpp -o CMakeFiles/2dSolverNS.dir/src/time_advance_RK3.cpp.s

CMakeFiles/2dSolverNS.dir/tests/computeErrorRHS.cpp.o: CMakeFiles/2dSolverNS.dir/flags.make
CMakeFiles/2dSolverNS.dir/tests/computeErrorRHS.cpp.o: ../tests/computeErrorRHS.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/vzendejasl/Documents/CaltechCourses/CS179/CS179Project/232aCFDCode_CUDA/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_10) "Building CXX object CMakeFiles/2dSolverNS.dir/tests/computeErrorRHS.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/2dSolverNS.dir/tests/computeErrorRHS.cpp.o -c /home/vzendejasl/Documents/CaltechCourses/CS179/CS179Project/232aCFDCode_CUDA/tests/computeErrorRHS.cpp

CMakeFiles/2dSolverNS.dir/tests/computeErrorRHS.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/2dSolverNS.dir/tests/computeErrorRHS.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/vzendejasl/Documents/CaltechCourses/CS179/CS179Project/232aCFDCode_CUDA/tests/computeErrorRHS.cpp > CMakeFiles/2dSolverNS.dir/tests/computeErrorRHS.cpp.i

CMakeFiles/2dSolverNS.dir/tests/computeErrorRHS.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/2dSolverNS.dir/tests/computeErrorRHS.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/vzendejasl/Documents/CaltechCourses/CS179/CS179Project/232aCFDCode_CUDA/tests/computeErrorRHS.cpp -o CMakeFiles/2dSolverNS.dir/tests/computeErrorRHS.cpp.s

# Object files for target 2dSolverNS
2dSolverNS_OBJECTS = \
"CMakeFiles/2dSolverNS.dir/main.cpp.o" \
"CMakeFiles/2dSolverNS.dir/src/compute_RHS.cpp.o" \
"CMakeFiles/2dSolverNS.dir/src/compute_curl.cpp.o" \
"CMakeFiles/2dSolverNS.dir/src/compute_projection_step.cpp.o" \
"CMakeFiles/2dSolverNS.dir/src/compute_RHS.cu.o" \
"CMakeFiles/2dSolverNS.dir/src/compute_projection_step.cu.o" \
"CMakeFiles/2dSolverNS.dir/src/time_advance_RK3.cu.o" \
"CMakeFiles/2dSolverNS.dir/src/helper_functions.cpp.o" \
"CMakeFiles/2dSolverNS.dir/src/time_advance_RK3.cpp.o" \
"CMakeFiles/2dSolverNS.dir/tests/computeErrorRHS.cpp.o"

# External object files for target 2dSolverNS
2dSolverNS_EXTERNAL_OBJECTS =

2dSolverNS: CMakeFiles/2dSolverNS.dir/main.cpp.o
2dSolverNS: CMakeFiles/2dSolverNS.dir/src/compute_RHS.cpp.o
2dSolverNS: CMakeFiles/2dSolverNS.dir/src/compute_curl.cpp.o
2dSolverNS: CMakeFiles/2dSolverNS.dir/src/compute_projection_step.cpp.o
2dSolverNS: CMakeFiles/2dSolverNS.dir/src/compute_RHS.cu.o
2dSolverNS: CMakeFiles/2dSolverNS.dir/src/compute_projection_step.cu.o
2dSolverNS: CMakeFiles/2dSolverNS.dir/src/time_advance_RK3.cu.o
2dSolverNS: CMakeFiles/2dSolverNS.dir/src/helper_functions.cpp.o
2dSolverNS: CMakeFiles/2dSolverNS.dir/src/time_advance_RK3.cpp.o
2dSolverNS: CMakeFiles/2dSolverNS.dir/tests/computeErrorRHS.cpp.o
2dSolverNS: CMakeFiles/2dSolverNS.dir/build.make
2dSolverNS: /usr/lib/x86_64-linux-gnu/libcufft.so
2dSolverNS: /usr/lib/x86_64-linux-gnu/libcublas.so
2dSolverNS: CMakeFiles/2dSolverNS.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/vzendejasl/Documents/CaltechCourses/CS179/CS179Project/232aCFDCode_CUDA/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_11) "Linking CXX executable 2dSolverNS"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/2dSolverNS.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/2dSolverNS.dir/build: 2dSolverNS

.PHONY : CMakeFiles/2dSolverNS.dir/build

CMakeFiles/2dSolverNS.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/2dSolverNS.dir/cmake_clean.cmake
.PHONY : CMakeFiles/2dSolverNS.dir/clean

CMakeFiles/2dSolverNS.dir/depend:
	cd /home/vzendejasl/Documents/CaltechCourses/CS179/CS179Project/232aCFDCode_CUDA/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/vzendejasl/Documents/CaltechCourses/CS179/CS179Project/232aCFDCode_CUDA /home/vzendejasl/Documents/CaltechCourses/CS179/CS179Project/232aCFDCode_CUDA /home/vzendejasl/Documents/CaltechCourses/CS179/CS179Project/232aCFDCode_CUDA/build /home/vzendejasl/Documents/CaltechCourses/CS179/CS179Project/232aCFDCode_CUDA/build /home/vzendejasl/Documents/CaltechCourses/CS179/CS179Project/232aCFDCode_CUDA/build/CMakeFiles/2dSolverNS.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/2dSolverNS.dir/depend
