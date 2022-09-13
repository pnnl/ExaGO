#!/bin/bash

# assuming running in exago root directory
spack install -j 4

spack module tcl refresh -y --delete-tree

spack module tcl loads -x exago &> ./buildsystem/clang-hip/crusher/dependencies.sh 
spack module tcl loads exago &> ./buildsystem/clang-hip/crusher/exago.sh 
chmod -R ugo+wrx /gpfs/alpine/proj-shared/csc359/cameron/spack-install/modules/cray-sles15-zen3
