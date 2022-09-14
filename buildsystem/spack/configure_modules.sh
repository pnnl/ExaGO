#!/bin/bash

# assuming running in exago root directory
# This assumes that we are installing from a binary mirror, and don't want to fetch files
spack install -j $1 && \

# This deletes the previous modules that are installed
# Either use a different module path than other users, or deal with frequent updates
# To use a different module path, you must update path in:
# 	- ./buildsystem/clang-hip/crusher/spack.yaml:25
#	- ./buildsystem/clang-hip/crusherVariables.sh
#	- Final chmod command at the end of this file
# TODO: - Test using environment variable
# TODO: - Hard code first line of file to be the module use /path/to/module
#  ^ This will enable the module path location to be specified here, and updated where necessary

# export SPACK_MODULE_INSTALL
# spack.yaml:
# module_root: $SPACK_MODULE_INSTALL
# crusherVariables.sh : remove module path
# dependencies.sh: add hard coded first line to include path

spack module tcl refresh -y --delete-tree && \
	
spack module tcl loads -r -x exago exago~full_optimizations &> ./buildsystem/clang-hip/crusher/dependencies.sh && \
spack module tcl loads -r -x exago exago+full_optimizations &> ./buildsystem/clang-hip/crusher/optimized-dependencies.sh && \
spack module tcl loads exago+full_optimizations &> ./buildsystem/clang-hip/crusher/exago-optimized.sh && \
spack module tcl loads exago~full_optimizations &> ./buildsystem/clang-hip/crusher/exago.sh && \

chmod -R ugo+wrx /gpfs/alpine/proj-shared/csc359/cameron/spack-install/test-modules/cray-sles15-zen3
