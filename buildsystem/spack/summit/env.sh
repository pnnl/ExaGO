#!/bin/bash

# Configure python
module load python/3.8-anaconda3

export SPACK_INSTALL=/gpfs/alpine/proj-shared/csc359/cameron/spack-install
export SPACK_MODULES=summit-modules
export SPACK_CACHE=/gpfs/alpine/proj-shared/csc359/$(whoami)/spack-cache
export SPACK_PYTHON=$OLCF_PYTHON_ANACONDA3_ROOT
export SPACK_DISABLE_LOCAL_CONFIG=1

