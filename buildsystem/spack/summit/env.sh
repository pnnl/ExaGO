#!/bin/bash

# Configure python
module load python/3.8-anaconda3

export SPACK_INSTALL=/gpfs/alpine2/stf006/world-shared/nkouk/parco-cuda11.4.2-spack-install
export SPACK_MODULES=summit-modules
export SPACK_CACHE=/gpfs/alpine2/stf006/world-shared/nkouk/parco-cuda11.4.2-spack-cache
export SPACK_PYTHON=$OLCF_PYTHON_ANACONDA3_ROOT
export SPACK_DISABLE_LOCAL_CONFIG=1

