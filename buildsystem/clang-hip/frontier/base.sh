#!/bin/bash

export MY_CLUSTER=frontier
export PROJ_DIR=/autofs/nccs-svm1_proj/eng151

module reset

# System modules
module load PrgEnv-gnu
module load craype-x86-trento
module load craype-accel-amd-gfx90a
module load rocm/6.3.1
module load cray-mpich
module load libfabric
module load cmake
module load cray-python

# The CC, CXX, FC environment variables can be overridden by Spack's module generation.
# The Frontier configuration was setup to avoid this override, but this is a
# potential pitfall (e.g., when configuring another system). 
# If all else fails, CMAKE_CXX_COMPILER supersedes these environment variables.
# Consider updating the build system to selectively mark HIP code in CMake.
export CC=/opt/rocm-6.3.1/llvm/bin/amdclang
export CXX=/opt/rocm-6.3.1/llvm/bin/amdclang++
export FC=/opt/rocm-6.3.1/llvm/bin/amdflang

export EXTRA_CMAKE_ARGS="$EXTRA_CMAKE_ARGS -DEXAGO_CTEST_LAUNCH_COMMAND='srun'"
export EXTRA_CMAKE_ARGS="$EXTRA_CMAKE_ARGS -DAMDGPU_TARGETS='gfx90a'"
