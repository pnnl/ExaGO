#!/bin/bash

module reset

# Configure python
module load cray-python/3.9.12.1 

# Configure base modules
module load PrgEnv-cray/8.3.3
module load craype-x86-trento
module load craype-accel-amd-gfx90a
module load cce/14.0.2
module load rocm/5.2.0
module load cray-mpich/8.1.23
module load libfabric/1.15.0.0

BASE=/ccs/proj/csc359/$(whoami)

export SPACK_INSTALL=/ccs/proj/csc359/spack-install
export SPACK_MODULES=modules
export SPACK_CACHE=$BASE/spack-cache
export SPACK_MIRROR=$BASE/spack-mirror
export SPACK_PYTHON=$CRAY_PYTHON_PREFIX
export SPACK_DISABLE_LOCAL_CONFIG=1

