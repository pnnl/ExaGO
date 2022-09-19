#!/bin/bash

SPACK_MIRROR=/gpfs/alpine/proj-shared/csc359/spack-install/source-mirror

spack mirror create --all --directory $SPACK_MIRROR
chmod -R ugo+rs $SPACK_MIRROR

