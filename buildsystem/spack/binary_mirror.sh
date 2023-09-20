#!/bin/bash

SPACK_MIRROR="${SPACK_MIRROR:?SPACK_MIRROR is unset. Please use the load_spack script first.}"

spack mirror create --all --directory $SPACK_MIRROR
chmod -R ugo+wrx $SPACK_MIRROR
