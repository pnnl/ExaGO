#!/bin/bash

SPACK_INSTALL="${SPACK_INSTALL:?SPACK_INSTALL is unset. Please use the load_spack script first.}"
SPACK_MIRROR=$SPACK_INSTALL/$SOURCE_CACHE

spack mirror create --all --directory $SPACK_MIRROR
chmod -R ugo+wrx $SPACK_MIRROR

