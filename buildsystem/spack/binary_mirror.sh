#!/bin/bash

SPACK_INSTALL="${SPACK_INSTALL:?SPACK_INSTALL is unset. Please use the load_spack script first.}"
SPACK_MIRROR="${SPACK_MIRROR:?SPACK_MIRROR is unset. Please use the load_spack script first.}"

echo "Spack source mirror is being created for you in $SPACK_MIRROR"

spack mirror create --all --directory $SPACK_MIRROR

