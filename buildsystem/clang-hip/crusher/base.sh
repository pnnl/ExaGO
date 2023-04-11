#!/bin/bash

export MY_CLUSTER=crusher
export PROJ_DIR=/autofs/nccs-svm1_proj/csc359

module reset

# System modules
module load PrgEnv-gnu-amd/8.3.3
module load amd-mixed/5.2.0
module load craype-x86-trento
module load craype-accel-amd-gfx90a
module load gcc/12.2.0
# Loaded through spack based module:
# - See ./buildsystem/spack/crusher/modules/dependencies.sh
#
# TODO: Figure out way to exclude MPI based on an abstract spec, as cray-mpich != openmpi
#	That way we could load MPI here. Also need to find way to only store this in one place...
#		- Perhaps through spack external find...
#		- Or spack scripting
#
# module load cray-mpich/8.1.23
module load libfabric/1.15.2.0

# Consider changing to $(which clang) as for deception
# While our stack is using the gcc toolchain,
# Using amdclang and the platform compilers is much smoother
# when invoking raw CMake. ABIs seem compatible
export CC=$(which amdclang)
export CXX=$(which amdclang++) 
export FC=$(which gfortran)

export EXTRA_CMAKE_ARGS="$EXTRA_CMAKE_ARGS -DEXAGO_CTEST_LAUNCH_COMMAND='srun'"

