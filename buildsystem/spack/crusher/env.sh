#!/bin/bash

module reset

# Configure python
module load cray-python/3.9.12.1 

# Configure base modules
module load PrgEnv-gnu/8.3.3
module load craype-x86-trento
module load craype-accel-amd-gfx90a
module load gcc/12.2.0
module load rocm/5.2.0
module load cray-mpich/8.1.23
module load libfabric/1.15.2.0

# Different scratch file systems depending on the cluster
# This is a simple fix, but we might have to re-factor if
# Frontier and Crusher end up being different enough
if sinfo | grep -q 'crusher'; then
  BASE=/gpfs/alpine/CSC359/scratch/$(whoami)
elif sinfo | grep -q 'frontier'; then
  BASE=/lustre/orion/csc359/scratch/$(whoami)
else
  echo "No Scratch filesystem available" && exit 1
fi

mkdir -p $BASE

# This is where the binaries and modules are installed
export SPACK_INSTALL=/ccs/proj/csc359/spack-install
# This is the prefix for the folder where modules are stored
export SPACK_MODULES=modules
export SPACK_CACHE=$BASE/spack-cache
export SPACK_MIRROR=$BASE/spack-mirror
export SPACK_PYTHON=$CRAY_PYTHON_PREFIX
export SPACK_DISABLE_LOCAL_CONFIG=1

