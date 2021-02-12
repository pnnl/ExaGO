#!/bin/bash

usage() {
  echo "Usage: ./build.sh [options]

--------------------------------------------------------------------------------

Long Description:

  This script is the entry point for ExaGO's continuous integration and default
  build configuration. The --build-only and --test-only options below build and
  test ExaGO with every option enabled. If you would like to build a smaller
  configuration, you will have to create a build directory and use the usual
  cmake workflow (eg edit variables in ccmake or pass command line arguments to
  cmake).

Clusters:

  By default, this script will attempt to determine the cluster it is being ran 
  on using the hostname command. If a known cluster is found, it's respective 
  script in the directory ./scripts/buildsystem will be sourced and the 
  variable MY_CLUSTER will be set. For example, on PNNL cluster Marianas, 
  hostname marianas.pnl.gov will be matched and 
  ./scripts/buildsystem/marianasVariables.sh will be sourced. If you would like 
  to add a cluster, create a script
  ./scripts/buildsystem/<my cluster>Variables.sh and specify the relevant
  environment variables. If the hostname is not correctly finding your cluster,
  you may specify MY_CLUSTER environment variable before running this script
  and the script will respect the environment variable. For example, on ORNL
  Ascent cluster, the hostname does not find the cluster, so we must specify
  MY_CLUSTER when running:

    $ MY_CLUSTER=ascent ./build.sh --build-only

Spack:

  Each supported variables script in ./scripts/buildsystem activates a spack
  environment with all dependencies configured. If you have built dependencies
  for ExaGO in a spack environment, you may simply activate the environment
  and run the build script specifying that you don't want to source any
  variables scripts, eg:

    $ MY_CLUSTER=none ./build.sh

--------------------------------------------------------------------------------

Options:

  --build-only    Only run the build stage of the script. This is useful for
                  local development.

  --test-only     Only run the test stage of the script. This should be ran
                  before every push to the repository or pull/merge request.
                  This run takes a significant amound of time.

--------------------------------------------------------------------------------

"
}

# This will be the catch-all trap handler after arguments are parsed.
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

makeArgs=${makeArgs:-"-j 8"}
ctestArgs=${ctestArgs:-"-VV"}
extraCmakeArgs=${extraCmakeArgs:-""}
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
  --help|-h)
    usage
    exit 0
    ;;
  *)
    echo "Argument $1 not recognized."
    usage
    exit 1
    ;;
  esac
done

set -xv
trap 'cleanup $? $LINENO' EXIT

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
  # We don't want all the shell functions we bring into scope to be printed out
  set -x
  source "scripts/buildsystem/$(echo $MY_CLUSTER)Variables.sh"
  set +x
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
