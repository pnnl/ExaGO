#!/bin/bash

# Define environment variables for where spack stores key files
# For now, SPACK_INSTALL is the path where everything spack related is installed
# If you want to modify the module install path, edit the spack.yaml manually
BASE=$(pwd)/spack_install
export SPACK_INSTALL=$BASE/binaries
export SPACK_MODULES=modules
export SPACK_CACHE=$(pwd)/spack-cache
export SPACK_DISABLE_LOCAL_CONFIG=1
export SPACK_MIRROR=$BASE/mirror
export SPACK_PYTHON=$(which python)

# Temporary directories used at build time
export tempdir=$SPACK_CACHE
export TMP=$SPACK_CACHE
export TMPDIR=$SPACK_CACHE

