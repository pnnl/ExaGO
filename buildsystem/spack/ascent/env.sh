#!/bin/bash

# Configure python
module load python/3.9-anaconda3

export SPACK_INSTALL=/gpfs/wolf/proj-shared/csc359/spack-install
export SPACK_MODULES=ascent-modules
export SPACK_CACHE=/gpfs/wolf/proj-shared/csc359/spack-cache/$(whoami)
export SPACK_PYTHON=$OLCF_PYTHON_ANACONDA3_ROOT
export SPACK_DISABLE_LOCAL_CONFIG=1
