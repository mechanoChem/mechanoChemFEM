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
CMAKE_COMMAND = /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/Cellar/cmake/3.7.2/bin/cmake

# The command to remove a file.
RM = /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/Cellar/cmake/3.7.2/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/wzhenlin/Github/mechanoChemFEM/build

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/wzhenlin/Github/mechanoChemFEM/build

# Include any dependencies generated for this target.
include CMakeFiles/main.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/main.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/main.dir/flags.make

CMakeFiles/main.dir/main.cc.o: CMakeFiles/main.dir/flags.make
CMakeFiles/main.dir/main.cc.o: main.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/wzhenlin/Github/mechanoChemFEM/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/main.dir/main.cc.o"
	/Applications/deal.II-8.5-brew.app/Contents/Resources/brew/bin/mpicxx   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/main.dir/main.cc.o -c /Users/wzhenlin/Github/mechanoChemFEM/build/main.cc

CMakeFiles/main.dir/main.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/main.dir/main.cc.i"
	/Applications/deal.II-8.5-brew.app/Contents/Resources/brew/bin/mpicxx  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/wzhenlin/Github/mechanoChemFEM/build/main.cc > CMakeFiles/main.dir/main.cc.i

CMakeFiles/main.dir/main.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/main.dir/main.cc.s"
	/Applications/deal.II-8.5-brew.app/Contents/Resources/brew/bin/mpicxx  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/wzhenlin/Github/mechanoChemFEM/build/main.cc -o CMakeFiles/main.dir/main.cc.s

CMakeFiles/main.dir/main.cc.o.requires:

.PHONY : CMakeFiles/main.dir/main.cc.o.requires

CMakeFiles/main.dir/main.cc.o.provides: CMakeFiles/main.dir/main.cc.o.requires
	$(MAKE) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/main.cc.o.provides.build
.PHONY : CMakeFiles/main.dir/main.cc.o.provides

CMakeFiles/main.dir/main.cc.o.provides.build: CMakeFiles/main.dir/main.cc.o


# Object files for target main
main_OBJECTS = \
"CMakeFiles/main.dir/main.cc.o"

# External object files for target main
main_EXTERNAL_OBJECTS =

main: CMakeFiles/main.dir/main.cc.o
main: CMakeFiles/main.dir/build.make
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/lib/libdeal_II.g.8.5.0.dylib
main: dealMultiPhysics.a
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/lib/libdeal_II.8.5.0.dylib
main: /usr/lib/libbz2.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libpike-blackbox.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libtrilinoscouplings.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libpanzer-disc-fe.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libpanzer-dof-mgr.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libpanzer-core.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libpiro.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/librol.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/librythmos.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libmuelu-adapters.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libmuelu-interface.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libmuelu.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libmoertel.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/liblocathyra.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/liblocaepetra.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/liblocalapack.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libloca.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libnoxepetra.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libnoxlapack.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libnox.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libphalanx.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libintrepid2.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libintrepid.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libteko.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libfei_trilinos.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libfei_base.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libstratimikos.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libstratimikosbelos.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libstratimikosaztecoo.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libstratimikosamesos.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libstratimikosml.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libstratimikosifpack.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libifpack2-adapters.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libifpack2.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libanasazitpetra.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libModeLaplace.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libanasaziepetra.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libanasazi.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libkomplex.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libamesos2.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libshylu.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libbelostpetra.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libbelosepetra.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libbelos.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libml.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libifpack.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libpamgen_extras.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libpamgen.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libamesos.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libgaleri-xpetra.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libgaleri-epetra.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libaztecoo.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libdpliris.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libisorropia.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/liboptipack.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libxpetra-sup.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libxpetra.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libthyratpetra.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libthyraepetraext.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libthyraepetra.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libthyracore.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libepetraext.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libtpetraext.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libtpetrainout.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libtpetra.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libkokkostsqr.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libtpetrakernels.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libtpetraclassiclinalg.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libtpetraclassicnodeapi.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libtpetraclassic.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libtriutils.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libglobipack.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libshards.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libzoltan.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libepetra.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libsacado.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/librtop.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libteuchoskokkoscomm.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libteuchoskokkoscompat.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libteuchosremainder.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libteuchosnumerics.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libteuchoscomm.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libteuchosparameterlist.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libteuchoscore.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libkokkosalgorithms.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libkokkoscontainers.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libkokkoscore.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libtpi.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libgtest.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libcholmod.a
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libamd.a
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libcolamd.a
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libsuitesparseconfig.a
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/opt/mumps/lib/libdmumps.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/opt/mumps/lib/libmumps_common.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/opt/mumps/lib/libpord.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libscalapack.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libsuperlu_dist.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libHYPRE.a
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libsuitesparseconfig.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libtbb.dylib
main: /usr/lib/libz.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/opt/parmetis/lib/libparmetis.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/opt/metis/lib/libmetis.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libptscotch.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libptscotcherr.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libscotch.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libscotcherr.dylib
main: /usr/lib/liblapack.dylib
main: /usr/lib/libblas.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libyaml-cpp.dylib
main: /usr/lib/libdl.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libhwloc.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/Cellar/open-mpi/2.1.0/lib/libmpi.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libumfpack.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libcholmod.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libccolamd.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libcolamd.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libcamd.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libamd.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libmetis.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libarpack.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libboost_iostreams-mt.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libboost_serialization-mt.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libboost_system-mt.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libboost_thread-mt.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libboost_regex-mt.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libboost_chrono-mt.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libboost_date_time-mt.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libboost_atomic-mt.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libgsl.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libgslcblas.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libhdf5_hl.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libhdf5.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libmuparser.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libnetcdf_c++.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libnetcdf.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libTKBO.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libTKBool.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libTKBRep.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libTKernel.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libTKFeat.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libTKFillet.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libTKG2d.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libTKG3d.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libTKGeomAlgo.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libTKGeomBase.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libTKHLR.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libTKIGES.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libTKMath.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libTKMesh.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libTKOffset.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libTKPrim.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libTKShHealing.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libTKSTEP.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libTKSTEPAttr.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libTKSTEPBase.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libTKSTEP209.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libTKSTL.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libTKTopAlgo.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libTKXSBase.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libp4est.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libsc.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libslepc.dylib
main: /Applications/deal.II-8.5-brew.app/Contents/Resources/brew/lib/libpetsc.dylib
main: CMakeFiles/main.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/wzhenlin/Github/mechanoChemFEM/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable main"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/main.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/main.dir/build: main

.PHONY : CMakeFiles/main.dir/build

CMakeFiles/main.dir/requires: CMakeFiles/main.dir/main.cc.o.requires

.PHONY : CMakeFiles/main.dir/requires

CMakeFiles/main.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/main.dir/cmake_clean.cmake
.PHONY : CMakeFiles/main.dir/clean

CMakeFiles/main.dir/depend:
	cd /Users/wzhenlin/Github/mechanoChemFEM/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/wzhenlin/Github/mechanoChemFEM/build /Users/wzhenlin/Github/mechanoChemFEM/build /Users/wzhenlin/Github/mechanoChemFEM/build /Users/wzhenlin/Github/mechanoChemFEM/build /Users/wzhenlin/Github/mechanoChemFEM/build/CMakeFiles/main.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/main.dir/depend

