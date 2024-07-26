#!/bin/bash

export MY_CLUSTER=crusher
export PROJ_DIR=/autofs/nccs-svm1_proj/eng145

module reset

# System modules
module load PrgEnv-amd
module load cpe/23.12
module load craype-x86-trento
module load craype-accel-amd-gfx90a
module load amd/5.7.1
module load rocm/5.7.1
module load cray-mpich/8.1.28
module load libfabric

# Consider changing to $(which clang) as for deception
export CC=/opt/rocm-5.7.1/llvm/bin/amdclang
export CXX=/opt/rocm-5.7.1/llvm/bin/amdclang++
export FC=/opt/rocm-5.7.1/llvm/bin/amdflang

export EXTRA_CMAKE_ARGS="$EXTRA_CMAKE_ARGS -DEXAGO_CTEST_LAUNCH_COMMAND='srun'"
export EXTRA_CMAKE_ARGS="$EXTRA_CMAKE_ARGS -DAMDGPU_TARGETS='gfx90a'"
