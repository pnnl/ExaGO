#!/bin/bash
# assuming running in exago root directory

SPACK_INSTALL="${SPACK_INSTALL:?SPACK_INSTALL is unset. Make sure to source load_spack.sh first}"

# Assuming MY_CLUSTER is configured
# TODO - fix this to be POSIX standardized with tolower across all scripts
base="./buildsystem/spack/${MY_CLUSTER,,""""}"

# Printing out loaded modules for debugging...
module list

# This assumes that we are installing from a binary mirror, and don't want to fetch files
# Make sure to use binary_mirror.sh if this is failing
spack develop --path=$(pwd) exago@develop && \
spack concretize -f && \
spack install --fail-fast -j $1 && \

# This deletes the previous modules that are installed
# Either use a different module path than other users, or deal with frequent updates
# To use a different module path, you must update the variable $SPACK_MODULES
# This is configured for the default for the branch in /buildsystem/spack/platform/env.sh
spack module tcl refresh -y && \
	
# We will create a new modules file, with the first line being the module path
# Note we redirect and destroy old file
# Only installing non-optimized version while config is not in spack develop
echo module use -a $SPACK_INSTALL/$SPACK_MODULES/$(spack arch) &> $base/modules/dependencies.sh && \
echo module use -a $SPACK_INSTALL/$SPACK_MODULES/$(spack arch) &> $base/modules/optimized-dependencies.sh && \
echo module use -a $SPACK_INSTALL/$SPACK_MODULES/$(spack arch) &> $base/modules/exago-optimized.sh && \
echo module use -a $SPACK_INSTALL/$SPACK_MODULES/$(spack arch) &> $base/modules/exago.sh && \

# Now we can append to the files
# It's a minor miracle that installing exago w/ different
# HiOp versions that different exago dev_path installs have
# different modules...
# TODO: would be to automatically get spack to pull the same
# 	git version as the branch being developed locally
spack module tcl loads -r -x exago build_type=RelWithDebInfo -x openssl exago &>> $base/modules/dependencies.sh && \
spack module tcl loads -r -x exago -x openssl exago build_type=Release &>> $base/modules/optimized-dependencies.sh && \
spack module tcl loads exago build_type=Release &>> $base/modules/exago-optimized.sh && \
spack module tcl loads exago build_type=RelWithDebInfo &>> $base/modules/exago.sh

exit_code=$?

# This makes the module and installation location globally readable which isn't ideal.
# Sticking to this avoids permission issues for other group members, but you
# should make sure that the install location is in a folder that is only
# accessible to group members
# Since we use this in CI and other users will create files we cannot chmod,
# We need to allow this command to fail
chmod -R ugo+wrx $SPACK_INSTALL > /dev/null 2>&1
chmod -R ugo+wrx $SPACK_CACHE > /dev/null 2>&1

# Should still change permissions before exiting, but should also return the exit
# code of spack related code
exit $exit_code
