#!/bin/bash

export MY_CLUSTER=crusher
export PROJ_DIR=/autofs/nccs-svm1_proj/csc359

module reset

# System modules
module load PrgEnv-amd
module load craype-accel-amd-gfx90a
module load amd/5.2.0
module load libfabric/1.15.0.0

# Consider changing to $(which clang) as for deception
export CC=/opt/rocm-5.2.0/llvm/bin/amdclang
export CXX=/opt/rocm-5.2.0/llvm/bin/amdclang++
export FC=/opt/rocm-5.2.0/llvm/bin/amdflang

export EXTRA_CMAKE_ARGS="$EXTRA_CMAKE_ARGS -DEXAGO_CTEST_LAUNCH_COMMAND='srun'"
