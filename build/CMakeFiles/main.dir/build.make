# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.8

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
CMAKE_COMMAND = /cvmfs/soft.computecanada.ca/nix/store/678hyhmplmrgwg96yyyxdvbpchb8maya-cmake-cursesUI-qt4UI-3.8.2/bin/cmake

# The command to remove a file.
RM = /cvmfs/soft.computecanada.ca/nix/store/678hyhmplmrgwg96yyyxdvbpchb8maya-cmake-cursesUI-qt4UI-3.8.2/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/aameri2/projects/rrg-nadaraja-ac/aameri2/Burgers-Equation

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/aameri2/projects/rrg-nadaraja-ac/aameri2/Burgers-Equation/build

# Include any dependencies generated for this target.
include CMakeFiles/main.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/main.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/main.dir/flags.make

CMakeFiles/main.dir/main.cc.o: CMakeFiles/main.dir/flags.make
CMakeFiles/main.dir/main.cc.o: ../main.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/aameri2/projects/rrg-nadaraja-ac/aameri2/Burgers-Equation/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/main.dir/main.cc.o"
	/cvmfs/soft.computecanada.ca/easybuild/software/2017/avx512/Compiler/intel2018.3/openmpi/3.1.2/bin/mpiCC  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/main.dir/main.cc.o -c /home/aameri2/projects/rrg-nadaraja-ac/aameri2/Burgers-Equation/main.cc

CMakeFiles/main.dir/main.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/main.dir/main.cc.i"
	/cvmfs/soft.computecanada.ca/easybuild/software/2017/avx512/Compiler/intel2018.3/openmpi/3.1.2/bin/mpiCC $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/aameri2/projects/rrg-nadaraja-ac/aameri2/Burgers-Equation/main.cc > CMakeFiles/main.dir/main.cc.i

CMakeFiles/main.dir/main.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/main.dir/main.cc.s"
	/cvmfs/soft.computecanada.ca/easybuild/software/2017/avx512/Compiler/intel2018.3/openmpi/3.1.2/bin/mpiCC $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/aameri2/projects/rrg-nadaraja-ac/aameri2/Burgers-Equation/main.cc -o CMakeFiles/main.dir/main.cc.s

CMakeFiles/main.dir/main.cc.o.requires:

.PHONY : CMakeFiles/main.dir/main.cc.o.requires

CMakeFiles/main.dir/main.cc.o.provides: CMakeFiles/main.dir/main.cc.o.requires
	$(MAKE) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/main.cc.o.provides.build
.PHONY : CMakeFiles/main.dir/main.cc.o.provides

CMakeFiles/main.dir/main.cc.o.provides.build: CMakeFiles/main.dir/main.cc.o


CMakeFiles/main.dir/util.cc.o: CMakeFiles/main.dir/flags.make
CMakeFiles/main.dir/util.cc.o: ../util.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/aameri2/projects/rrg-nadaraja-ac/aameri2/Burgers-Equation/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/main.dir/util.cc.o"
	/cvmfs/soft.computecanada.ca/easybuild/software/2017/avx512/Compiler/intel2018.3/openmpi/3.1.2/bin/mpiCC  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/main.dir/util.cc.o -c /home/aameri2/projects/rrg-nadaraja-ac/aameri2/Burgers-Equation/util.cc

CMakeFiles/main.dir/util.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/main.dir/util.cc.i"
	/cvmfs/soft.computecanada.ca/easybuild/software/2017/avx512/Compiler/intel2018.3/openmpi/3.1.2/bin/mpiCC $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/aameri2/projects/rrg-nadaraja-ac/aameri2/Burgers-Equation/util.cc > CMakeFiles/main.dir/util.cc.i

CMakeFiles/main.dir/util.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/main.dir/util.cc.s"
	/cvmfs/soft.computecanada.ca/easybuild/software/2017/avx512/Compiler/intel2018.3/openmpi/3.1.2/bin/mpiCC $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/aameri2/projects/rrg-nadaraja-ac/aameri2/Burgers-Equation/util.cc -o CMakeFiles/main.dir/util.cc.s

CMakeFiles/main.dir/util.cc.o.requires:

.PHONY : CMakeFiles/main.dir/util.cc.o.requires

CMakeFiles/main.dir/util.cc.o.provides: CMakeFiles/main.dir/util.cc.o.requires
	$(MAKE) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/util.cc.o.provides.build
.PHONY : CMakeFiles/main.dir/util.cc.o.provides

CMakeFiles/main.dir/util.cc.o.provides.build: CMakeFiles/main.dir/util.cc.o


CMakeFiles/main.dir/parameters.cc.o: CMakeFiles/main.dir/flags.make
CMakeFiles/main.dir/parameters.cc.o: ../parameters.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/aameri2/projects/rrg-nadaraja-ac/aameri2/Burgers-Equation/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/main.dir/parameters.cc.o"
	/cvmfs/soft.computecanada.ca/easybuild/software/2017/avx512/Compiler/intel2018.3/openmpi/3.1.2/bin/mpiCC  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/main.dir/parameters.cc.o -c /home/aameri2/projects/rrg-nadaraja-ac/aameri2/Burgers-Equation/parameters.cc

CMakeFiles/main.dir/parameters.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/main.dir/parameters.cc.i"
	/cvmfs/soft.computecanada.ca/easybuild/software/2017/avx512/Compiler/intel2018.3/openmpi/3.1.2/bin/mpiCC $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/aameri2/projects/rrg-nadaraja-ac/aameri2/Burgers-Equation/parameters.cc > CMakeFiles/main.dir/parameters.cc.i

CMakeFiles/main.dir/parameters.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/main.dir/parameters.cc.s"
	/cvmfs/soft.computecanada.ca/easybuild/software/2017/avx512/Compiler/intel2018.3/openmpi/3.1.2/bin/mpiCC $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/aameri2/projects/rrg-nadaraja-ac/aameri2/Burgers-Equation/parameters.cc -o CMakeFiles/main.dir/parameters.cc.s

CMakeFiles/main.dir/parameters.cc.o.requires:

.PHONY : CMakeFiles/main.dir/parameters.cc.o.requires

CMakeFiles/main.dir/parameters.cc.o.provides: CMakeFiles/main.dir/parameters.cc.o.requires
	$(MAKE) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/parameters.cc.o.provides.build
.PHONY : CMakeFiles/main.dir/parameters.cc.o.provides

CMakeFiles/main.dir/parameters.cc.o.provides.build: CMakeFiles/main.dir/parameters.cc.o


CMakeFiles/main.dir/equations.cc.o: CMakeFiles/main.dir/flags.make
CMakeFiles/main.dir/equations.cc.o: ../equations.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/aameri2/projects/rrg-nadaraja-ac/aameri2/Burgers-Equation/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/main.dir/equations.cc.o"
	/cvmfs/soft.computecanada.ca/easybuild/software/2017/avx512/Compiler/intel2018.3/openmpi/3.1.2/bin/mpiCC  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/main.dir/equations.cc.o -c /home/aameri2/projects/rrg-nadaraja-ac/aameri2/Burgers-Equation/equations.cc

CMakeFiles/main.dir/equations.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/main.dir/equations.cc.i"
	/cvmfs/soft.computecanada.ca/easybuild/software/2017/avx512/Compiler/intel2018.3/openmpi/3.1.2/bin/mpiCC $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/aameri2/projects/rrg-nadaraja-ac/aameri2/Burgers-Equation/equations.cc > CMakeFiles/main.dir/equations.cc.i

CMakeFiles/main.dir/equations.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/main.dir/equations.cc.s"
	/cvmfs/soft.computecanada.ca/easybuild/software/2017/avx512/Compiler/intel2018.3/openmpi/3.1.2/bin/mpiCC $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/aameri2/projects/rrg-nadaraja-ac/aameri2/Burgers-Equation/equations.cc -o CMakeFiles/main.dir/equations.cc.s

CMakeFiles/main.dir/equations.cc.o.requires:

.PHONY : CMakeFiles/main.dir/equations.cc.o.requires

CMakeFiles/main.dir/equations.cc.o.provides: CMakeFiles/main.dir/equations.cc.o.requires
	$(MAKE) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/equations.cc.o.provides.build
.PHONY : CMakeFiles/main.dir/equations.cc.o.provides

CMakeFiles/main.dir/equations.cc.o.provides.build: CMakeFiles/main.dir/equations.cc.o


CMakeFiles/main.dir/problem.cc.o: CMakeFiles/main.dir/flags.make
CMakeFiles/main.dir/problem.cc.o: ../problem.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/aameri2/projects/rrg-nadaraja-ac/aameri2/Burgers-Equation/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMakeFiles/main.dir/problem.cc.o"
	/cvmfs/soft.computecanada.ca/easybuild/software/2017/avx512/Compiler/intel2018.3/openmpi/3.1.2/bin/mpiCC  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/main.dir/problem.cc.o -c /home/aameri2/projects/rrg-nadaraja-ac/aameri2/Burgers-Equation/problem.cc

CMakeFiles/main.dir/problem.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/main.dir/problem.cc.i"
	/cvmfs/soft.computecanada.ca/easybuild/software/2017/avx512/Compiler/intel2018.3/openmpi/3.1.2/bin/mpiCC $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/aameri2/projects/rrg-nadaraja-ac/aameri2/Burgers-Equation/problem.cc > CMakeFiles/main.dir/problem.cc.i

CMakeFiles/main.dir/problem.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/main.dir/problem.cc.s"
	/cvmfs/soft.computecanada.ca/easybuild/software/2017/avx512/Compiler/intel2018.3/openmpi/3.1.2/bin/mpiCC $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/aameri2/projects/rrg-nadaraja-ac/aameri2/Burgers-Equation/problem.cc -o CMakeFiles/main.dir/problem.cc.s

CMakeFiles/main.dir/problem.cc.o.requires:

.PHONY : CMakeFiles/main.dir/problem.cc.o.requires

CMakeFiles/main.dir/problem.cc.o.provides: CMakeFiles/main.dir/problem.cc.o.requires
	$(MAKE) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/problem.cc.o.provides.build
.PHONY : CMakeFiles/main.dir/problem.cc.o.provides

CMakeFiles/main.dir/problem.cc.o.provides.build: CMakeFiles/main.dir/problem.cc.o


# Object files for target main
main_OBJECTS = \
"CMakeFiles/main.dir/main.cc.o" \
"CMakeFiles/main.dir/util.cc.o" \
"CMakeFiles/main.dir/parameters.cc.o" \
"CMakeFiles/main.dir/equations.cc.o" \
"CMakeFiles/main.dir/problem.cc.o"

# External object files for target main
main_EXTERNAL_OBJECTS =

main: CMakeFiles/main.dir/main.cc.o
main: CMakeFiles/main.dir/util.cc.o
main: CMakeFiles/main.dir/parameters.cc.o
main: CMakeFiles/main.dir/equations.cc.o
main: CMakeFiles/main.dir/problem.cc.o
main: CMakeFiles/main.dir/build.make
main: /cvmfs/soft.computecanada.ca/easybuild/software/2017/avx512/MPI/intel2018.3/openmpi3.1/dealii/8.5.0/lib/libdeal_II.so.8.5.0
main: /cvmfs/soft.computecanada.ca/nix/var/nix/profiles/16.09/lib/libbz2.so
main: /cvmfs/soft.computecanada.ca/easybuild/software/2017/avx512/Compiler/intel2018.3/openmpi/3.1.2/lib/libmpi_usempif08.so
main: /cvmfs/soft.computecanada.ca/easybuild/software/2017/avx512/Compiler/intel2018.3/openmpi/3.1.2/lib/libmpi_usempi_ignore_tkr.so
main: /cvmfs/soft.computecanada.ca/easybuild/software/2017/avx512/Compiler/intel2018.3/openmpi/3.1.2/lib/libmpi_mpifh.so
main: /cvmfs/soft.computecanada.ca/nix/var/nix/profiles/16.09/lib/libz.so
main: /cvmfs/soft.computecanada.ca/easybuild/software/2017/avx512/MPI/intel2018.3/openmpi3.1/trilinos/12.10.1/lib/libmuelu-adapters.so
main: /cvmfs/soft.computecanada.ca/easybuild/software/2017/avx512/MPI/intel2018.3/openmpi3.1/trilinos/12.10.1/lib/libmuelu-interface.so
main: /cvmfs/soft.computecanada.ca/easybuild/software/2017/avx512/MPI/intel2018.3/openmpi3.1/trilinos/12.10.1/lib/libmuelu.so
main: /cvmfs/soft.computecanada.ca/easybuild/software/2017/avx512/MPI/intel2018.3/openmpi3.1/trilinos/12.10.1/lib/libteko.so
main: /cvmfs/soft.computecanada.ca/easybuild/software/2017/avx512/MPI/intel2018.3/openmpi3.1/trilinos/12.10.1/lib/libstratimikos.so
main: /cvmfs/soft.computecanada.ca/easybuild/software/2017/avx512/MPI/intel2018.3/openmpi3.1/trilinos/12.10.1/lib/libstratimikosbelos.so
main: /cvmfs/soft.computecanada.ca/easybuild/software/2017/avx512/MPI/intel2018.3/openmpi3.1/trilinos/12.10.1/lib/libstratimikosaztecoo.so
main: /cvmfs/soft.computecanada.ca/easybuild/software/2017/avx512/MPI/intel2018.3/openmpi3.1/trilinos/12.10.1/lib/libstratimikosamesos.so
main: /cvmfs/soft.computecanada.ca/easybuild/software/2017/avx512/MPI/intel2018.3/openmpi3.1/trilinos/12.10.1/lib/libstratimikosml.so
main: /cvmfs/soft.computecanada.ca/easybuild/software/2017/avx512/MPI/intel2018.3/openmpi3.1/trilinos/12.10.1/lib/libstratimikosifpack.so
main: /cvmfs/soft.computecanada.ca/easybuild/software/2017/avx512/MPI/intel2018.3/openmpi3.1/trilinos/12.10.1/lib/libifpack2-adapters.so
main: /cvmfs/soft.computecanada.ca/easybuild/software/2017/avx512/MPI/intel2018.3/openmpi3.1/trilinos/12.10.1/lib/libifpack2.so
main: /cvmfs/soft.computecanada.ca/easybuild/software/2017/avx512/MPI/intel2018.3/openmpi3.1/trilinos/12.10.1/lib/libanasazitpetra.so
main: /cvmfs/soft.computecanada.ca/easybuild/software/2017/avx512/MPI/intel2018.3/openmpi3.1/trilinos/12.10.1/lib/libModeLaplace.so
main: /cvmfs/soft.computecanada.ca/easybuild/software/2017/avx512/MPI/intel2018.3/openmpi3.1/trilinos/12.10.1/lib/libanasaziepetra.so
main: /cvmfs/soft.computecanada.ca/easybuild/software/2017/avx512/MPI/intel2018.3/openmpi3.1/trilinos/12.10.1/lib/libanasazi.so
main: /cvmfs/soft.computecanada.ca/easybuild/software/2017/avx512/MPI/intel2018.3/openmpi3.1/trilinos/12.10.1/lib/libamesos2.so
main: /cvmfs/soft.computecanada.ca/easybuild/software/2017/avx512/MPI/intel2018.3/openmpi3.1/trilinos/12.10.1/lib/libbelostpetra.so
main: /cvmfs/soft.computecanada.ca/easybuild/software/2017/avx512/MPI/intel2018.3/openmpi3.1/trilinos/12.10.1/lib/libbelosepetra.so
main: /cvmfs/soft.computecanada.ca/easybuild/software/2017/avx512/MPI/intel2018.3/openmpi3.1/trilinos/12.10.1/lib/libbelos.so
main: /cvmfs/soft.computecanada.ca/easybuild/software/2017/avx512/MPI/intel2018.3/openmpi3.1/trilinos/12.10.1/lib/libml.so
main: /cvmfs/soft.computecanada.ca/easybuild/software/2017/avx512/MPI/intel2018.3/openmpi3.1/trilinos/12.10.1/lib/libifpack.so
main: /cvmfs/soft.computecanada.ca/easybuild/software/2017/avx512/MPI/intel2018.3/openmpi3.1/trilinos/12.10.1/lib/libzoltan2.so
main: /cvmfs/soft.computecanada.ca/easybuild/software/2017/avx512/MPI/intel2018.3/openmpi3.1/trilinos/12.10.1/lib/libpamgen_extras.so
main: /cvmfs/soft.computecanada.ca/easybuild/software/2017/avx512/MPI/intel2018.3/openmpi3.1/trilinos/12.10.1/lib/libpamgen.so
main: /cvmfs/soft.computecanada.ca/easybuild/software/2017/avx512/MPI/intel2018.3/openmpi3.1/trilinos/12.10.1/lib/libamesos.so
main: /cvmfs/soft.computecanada.ca/easybuild/software/2017/avx512/MPI/intel2018.3/openmpi3.1/trilinos/12.10.1/lib/libgaleri-xpetra.so
main: /cvmfs/soft.computecanada.ca/easybuild/software/2017/avx512/MPI/intel2018.3/openmpi3.1/trilinos/12.10.1/lib/libgaleri-epetra.so
main: /cvmfs/soft.computecanada.ca/easybuild/software/2017/avx512/MPI/intel2018.3/openmpi3.1/trilinos/12.10.1/lib/libaztecoo.so
main: /cvmfs/soft.computecanada.ca/easybuild/software/2017/avx512/MPI/intel2018.3/openmpi3.1/trilinos/12.10.1/lib/libisorropia.so
main: /cvmfs/soft.computecanada.ca/easybuild/software/2017/avx512/MPI/intel2018.3/openmpi3.1/trilinos/12.10.1/lib/libxpetra-sup.so
main: /cvmfs/soft.computecanada.ca/easybuild/software/2017/avx512/MPI/intel2018.3/openmpi3.1/trilinos/12.10.1/lib/libxpetra.so
main: /cvmfs/soft.computecanada.ca/easybuild/software/2017/avx512/MPI/intel2018.3/openmpi3.1/trilinos/12.10.1/lib/libthyratpetra.so
main: /cvmfs/soft.computecanada.ca/easybuild/software/2017/avx512/MPI/intel2018.3/openmpi3.1/trilinos/12.10.1/lib/libthyraepetraext.so
main: /cvmfs/soft.computecanada.ca/easybuild/software/2017/avx512/MPI/intel2018.3/openmpi3.1/trilinos/12.10.1/lib/libthyraepetra.so
main: /cvmfs/soft.computecanada.ca/easybuild/software/2017/avx512/MPI/intel2018.3/openmpi3.1/trilinos/12.10.1/lib/libthyracore.so
main: /cvmfs/soft.computecanada.ca/easybuild/software/2017/avx512/MPI/intel2018.3/openmpi3.1/trilinos/12.10.1/lib/libepetraext.so
main: /cvmfs/soft.computecanada.ca/easybuild/software/2017/avx512/MPI/intel2018.3/openmpi3.1/trilinos/12.10.1/lib/libtpetraext.so
main: /cvmfs/soft.computecanada.ca/easybuild/software/2017/avx512/MPI/intel2018.3/openmpi3.1/trilinos/12.10.1/lib/libtpetrainout.so
main: /cvmfs/soft.computecanada.ca/easybuild/software/2017/avx512/MPI/intel2018.3/openmpi3.1/trilinos/12.10.1/lib/libtpetra.so
main: /cvmfs/soft.computecanada.ca/easybuild/software/2017/avx512/MPI/intel2018.3/openmpi3.1/trilinos/12.10.1/lib/libkokkostsqr.so
main: /cvmfs/soft.computecanada.ca/easybuild/software/2017/avx512/MPI/intel2018.3/openmpi3.1/trilinos/12.10.1/lib/libtpetrakernels.so
main: /cvmfs/soft.computecanada.ca/easybuild/software/2017/avx512/MPI/intel2018.3/openmpi3.1/trilinos/12.10.1/lib/libtpetraclassiclinalg.so
main: /cvmfs/soft.computecanada.ca/easybuild/software/2017/avx512/MPI/intel2018.3/openmpi3.1/trilinos/12.10.1/lib/libtpetraclassicnodeapi.so
main: /cvmfs/soft.computecanada.ca/easybuild/software/2017/avx512/MPI/intel2018.3/openmpi3.1/trilinos/12.10.1/lib/libtpetraclassic.so
main: /cvmfs/soft.computecanada.ca/easybuild/software/2017/avx512/MPI/intel2018.3/openmpi3.1/trilinos/12.10.1/lib/libtriutils.so
main: /cvmfs/soft.computecanada.ca/easybuild/software/2017/avx512/MPI/intel2018.3/openmpi3.1/trilinos/12.10.1/lib/libzoltan.so
main: /cvmfs/soft.computecanada.ca/easybuild/software/2017/avx512/MPI/intel2018.3/openmpi3.1/trilinos/12.10.1/lib/libepetra.so
main: /cvmfs/soft.computecanada.ca/easybuild/software/2017/avx512/MPI/intel2018.3/openmpi3.1/trilinos/12.10.1/lib/libsacado.so
main: /cvmfs/soft.computecanada.ca/easybuild/software/2017/avx512/MPI/intel2018.3/openmpi3.1/trilinos/12.10.1/lib/librtop.so
main: /cvmfs/soft.computecanada.ca/easybuild/software/2017/avx512/MPI/intel2018.3/openmpi3.1/trilinos/12.10.1/lib/libteuchoskokkoscomm.so
main: /cvmfs/soft.computecanada.ca/easybuild/software/2017/avx512/MPI/intel2018.3/openmpi3.1/trilinos/12.10.1/lib/libteuchoskokkoscompat.so
main: /cvmfs/soft.computecanada.ca/easybuild/software/2017/avx512/MPI/intel2018.3/openmpi3.1/trilinos/12.10.1/lib/libteuchosremainder.so
main: /cvmfs/soft.computecanada.ca/easybuild/software/2017/avx512/MPI/intel2018.3/openmpi3.1/trilinos/12.10.1/lib/libteuchosnumerics.so
main: /cvmfs/soft.computecanada.ca/easybuild/software/2017/avx512/MPI/intel2018.3/openmpi3.1/trilinos/12.10.1/lib/libteuchoscomm.so
main: /cvmfs/soft.computecanada.ca/easybuild/software/2017/avx512/MPI/intel2018.3/openmpi3.1/trilinos/12.10.1/lib/libteuchosparameterlist.so
main: /cvmfs/soft.computecanada.ca/easybuild/software/2017/avx512/MPI/intel2018.3/openmpi3.1/trilinos/12.10.1/lib/libteuchoscore.so
main: /cvmfs/soft.computecanada.ca/easybuild/software/2017/avx512/MPI/intel2018.3/openmpi3.1/trilinos/12.10.1/lib/libkokkosalgorithms.so
main: /cvmfs/soft.computecanada.ca/easybuild/software/2017/avx512/MPI/intel2018.3/openmpi3.1/trilinos/12.10.1/lib/libkokkoscontainers.so
main: /cvmfs/soft.computecanada.ca/easybuild/software/2017/avx512/MPI/intel2018.3/openmpi3.1/trilinos/12.10.1/lib/libkokkoscore.so
main: /cvmfs/soft.computecanada.ca/easybuild/software/2017/Core/imkl/2018.3.222/mkl/lib/intel64/libmkl_rt.so
main: /cvmfs/soft.computecanada.ca/easybuild/software/2017/avx512/Compiler/intel2018.3/openmpi/3.1.2/lib/libmpi_cxx.so
main: /cvmfs/soft.computecanada.ca/easybuild/software/2017/avx512/Compiler/intel2018.3/gsl/2.5/lib/libgsl.a
main: /cvmfs/soft.computecanada.ca/easybuild/software/2017/avx512/MPI/intel2018.3/openmpi3.1/p4est/1.1/lib/libp4est.so
main: /cvmfs/soft.computecanada.ca/easybuild/software/2017/avx512/MPI/intel2018.3/openmpi3.1/p4est/1.1/lib/libsc.so
main: /cvmfs/soft.computecanada.ca/easybuild/software/2017/avx512/Compiler/intel2018.3/openmpi/3.1.2/lib/libmpi.so
main: CMakeFiles/main.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/aameri2/projects/rrg-nadaraja-ac/aameri2/Burgers-Equation/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Linking CXX executable main"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/main.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/main.dir/build: main

.PHONY : CMakeFiles/main.dir/build

CMakeFiles/main.dir/requires: CMakeFiles/main.dir/main.cc.o.requires
CMakeFiles/main.dir/requires: CMakeFiles/main.dir/util.cc.o.requires
CMakeFiles/main.dir/requires: CMakeFiles/main.dir/parameters.cc.o.requires
CMakeFiles/main.dir/requires: CMakeFiles/main.dir/equations.cc.o.requires
CMakeFiles/main.dir/requires: CMakeFiles/main.dir/problem.cc.o.requires

.PHONY : CMakeFiles/main.dir/requires

CMakeFiles/main.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/main.dir/cmake_clean.cmake
.PHONY : CMakeFiles/main.dir/clean

CMakeFiles/main.dir/depend:
	cd /home/aameri2/projects/rrg-nadaraja-ac/aameri2/Burgers-Equation/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/aameri2/projects/rrg-nadaraja-ac/aameri2/Burgers-Equation /home/aameri2/projects/rrg-nadaraja-ac/aameri2/Burgers-Equation /home/aameri2/projects/rrg-nadaraja-ac/aameri2/Burgers-Equation/build /home/aameri2/projects/rrg-nadaraja-ac/aameri2/Burgers-Equation/build /home/aameri2/projects/rrg-nadaraja-ac/aameri2/Burgers-Equation/build/CMakeFiles/main.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/main.dir/depend

