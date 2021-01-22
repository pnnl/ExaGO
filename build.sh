#!/bin/bash

cleanup() {
  echo
  echo Exit code $1 caught in build script.
  echo
  if [[ "$1" == "0" ]]; then
    echo BUILD_STATUS:0
  else
    echo
    echo Failure found on line $2 in build script.
    echo
    echo BUILD_STATUS:1
  fi
}

trap 'cleanup $? $LINENO' EXIT

set -xv
makeArgs=${makeArgs:-"-j 8"}
ctestArgs=${ctestArgs:-"-VV"}
extraCmakeArgs=${extraCmakeArgs:-""}
export CXXFLAGS='-I/share/apps/magma/2.5.2/cuda10.2/include -I/share/apps/cuda/10.2/include/'
export OMPI_MCA_btl="^vader,tcp,openib,uct"
export BUILD=${BUILD:-1}
export TEST=${TEST:-1}
export srcdir=${srcdir:-$PWD}
export builddir=${builddir:-$(pwd)/build}
export installdir=${installdir:-$(pwd)/install}
export BUILD_MATRIX=${BUILD_MATRIX:-0}

while [[ $# -gt 0 ]]; do
  case $1 in
  --short)
    echo "Only running representative subset of unit tests."
    ctestArgs="$ctestArgs -R UNIT_TEST"
    shift
    ;;
  --build-only)
    export TEST=0 BUILD=1
    shift
    ;;
  --test-only)
    export TEST=1 BUILD=0
    shift
    ;;
  --matrix)
    export BUILD_MATRIX=1
    shift
    ;;
  --ci)
    # This option is only useful if you specify all indices for the build
    # matrix options, such as CI_RAJAOPT. See
    # $EXAGODIR/scripts/buildsystem/buildMatrix.sh for more details.
    export BUILD_MATRIX=1
    export BUILD_MATRIX_PARALLEL=1
    shift
    ;;
  *)
    echo "Argument $1 not recognized."
    exit 1
    ;;
  esac
done

if [[ ! -v MY_CLUSTER ]]
then
  export MY_CLUSTER=`uname -n | sed -e 's/[0-9]//g' -e 's/\..*//'`
fi

# Correctly identify clusters based on hostname
case $MY_CLUSTER in
  newell*)
    export MY_CLUSTER=newell
    ;;
  dl*|marianas|*fat*)
    export MY_CLUSTER=marianas
    ;;
  *)
    echo "Cluster $MY_CLUSTER not identified - you'll have to set relevant variables manually."
    ;;
esac

ulimit -s unlimited || echo 'Could not set stack size to unlimited.'
ulimit -l unlimited || echo 'Could not set max locked memory to unlimited.'

. /etc/profile.d/modules.sh
module purge

if [[ -f "scripts/buildsystem/$(echo $MY_CLUSTER)Variables.sh" ]]; then
  source "scripts/buildsystem/$(echo $MY_CLUSTER)Variables.sh"
  echo Sourced system-specific variables for $MY_CLUSTER
fi

export LD_LIBRARY_PATH="$magma_dir/lib:$hiop_dir/lib:$LD_LIBRARY_PATH"
export CTEST_OUTPUT_ON_FAILURE=1

module list

#  NOTE: The following is required when running from Gitlab CI via slurm job
base_path=`dirname $0`
if [ -z "$SLURM_SUBMIT_DIR" ]; then
  cd $base_path          || exit 1
fi

if [[ $BUILD_MATRIX -eq 1 ]]; then
  echo Running build matrix
  source scripts/buildsystem/buildMatrix.sh
  buildMatrix
  exit $?
else
  echo Running defualt build
  source scripts/buildsystem/defaultBuild.sh
  defaultBuild
  exit $?
fi
