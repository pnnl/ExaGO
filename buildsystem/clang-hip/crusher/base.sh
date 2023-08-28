#!/bin/bash

export MY_CLUSTER=crusher
export PROJ_DIR=/autofs/nccs-svm1_proj/csc359

module purge

# System modules
module load rocm/5.2.0
module load libfabric/1.15.0.0

# Consider changing to $(which clang) as for deception
export CC=/opt/rocm-5.2.0/llvm/bin/clang
export CXX=/opt/rocm-5.2.0/llvm/bin/clang++
export FC=/opt/rocm-5.2.0/llvm/bin/flang

export EXTRA_CMAKE_ARGS="$EXTRA_CMAKE_ARGS -DEXAGO_CTEST_LAUNCH_COMMAND='srun'"
export EXTRA_CMAKE_ARGS="$EXTRA_CMAKE_ARGS -DAMDGPU_TARGETS='gfx90a'"

