#!/bin/bash

# Configure python
module load cray-python/3.9.12.1 

# Configure base modules
module load PrgEnv-amd
module load craype-x86-trento
module load craype-accel-amd-gfx90a
module load amd/5.2.0
module load cray-mpich/8.1.23
module load libfabric

export SPACK_INSTALL=$PROJWORK/csc359/$(whoami)/spack-install
export SPACK_MODULES=test-modules
export SPACK_CACHE=$PROJWORK/csc359/$(whoami)/spack-cache
export SPACK_PYTHON=$CRAY_PYTHON_PREFIX
export SPACK_DISABLE_LOCAL_CONFIG=1

