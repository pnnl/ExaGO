#!/bin/bash

export MY_CLUSTER=frontier
export PROJ_DIR=/autofs/nccs-svm1_proj/eng151

module reset

# System modules
module load PrgEnv-amd
module load craype-x86-trento
module load craype-accel-amd-gfx90a
module load amd/6.3.1
module load rocm/6.3.1
module load cray-mpich
module load libfabric

# Consider changing to $(which clang) as for deception
export CC=/opt/rocm-6.3.1/llvm/bin/amdclang
export CXX=/opt/rocm-6.3.1/llvm/bin/amdclang++
export FC=/opt/rocm-6.3.1/llvm/bin/amdflang

export EXTRA_CMAKE_ARGS="$EXTRA_CMAKE_ARGS -DEXAGO_CTEST_LAUNCH_COMMAND='srun'"
export EXTRA_CMAKE_ARGS="$EXTRA_CMAKE_ARGS -DAMDGPU_TARGETS='gfx90a'"
