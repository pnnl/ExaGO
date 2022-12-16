#!/bin/bash

SPACK_MIRROR="${SPACK_MIRROR:?SPACK_MIRROR is unset. Please use the load_spack script first.}"

echo "Spack source mirror is being created for you in $SPACK_MIRROR"

spack mirror create --all --directory $SPACK_MIRROR

