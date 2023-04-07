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
module load cray-mpich/8.1.23
module load libfabric/1.15.2.0

# Consider changing to $(which clang) as for deception
# While our stack is using the gcc toolchain,
# Using amdclang and the platform compilers is much smoother
# when invoking raw CMake. ABIs seem compatible
export CC=$(which amdclang)
export CXX=$(which amdclang++) 
export FC=$(which gfortran)

export EXTRA_CMAKE_ARGS="$EXTRA_CMAKE_ARGS -DEXAGO_CTEST_LAUNCH_COMMAND='srun'"

