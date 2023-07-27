#!/bin/bash

export MY_CLUSTER=crusher
export PROJ_DIR=/autofs/nccs-svm1_proj/csc359

module reset

# System modules
module load PrgEnv-amd
module load craype-x86-trento
module load craype-accel-amd-gfx90a
module load amd/5.2.0
module load cray-mpich/8.1.25
module load libfabric

# Consider changing to $(which clang) as for deception
export CC=/opt/rocm-5.2.0/llvm/bin/amdclang
export CXX=/opt/rocm-5.2.0/llvm/bin/amdclang++
export FC=/opt/rocm-5.2.0/llvm/bin/amdflang

export EXTRA_CMAKE_ARGS="$EXTRA_CMAKE_ARGS -DEXAGO_CTEST_LAUNCH_COMMAND='srun'"
export EXTRA_CMAKE_ARGS="$EXTRA_CMAKE_ARGS -DAMDGPU_TARGETS='gfx90a'"

