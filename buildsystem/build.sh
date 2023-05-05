#!/bin/bash

exit() {
  # Clear all trap handlers so this isn't echo'ed multiple times, potentially
  # throwing off the CI script watching for this output
  trap - `seq 1 31`

  # If called without an argument, assume not an error
  local ec=${1:-0}

  # Echo the snippet the CI script is looking for
  echo BUILD_STATUS:${ec}

  # Actually exit with that code, although it won't matter in most cases, as CI
  # is only looking for the string 'BUILD_STATUS:N'
  builtin exit ${ec}
}

# This will be the catch-all trap handler after arguments are parsed.
cleanup() {
  # Clear all trap handlers
  trap - `seq 1 31`

  # When 'trap' is invoked, each signal handler will be a curried version of
  # this function which has the first argument bound to the signal it's catching
  local sig=$1

  echo
  echo Exit code $2 caught in build script triggered by signal ${sig}.
  echo

  exit $2
}

if [[ ! -f $PWD/buildsystem/build.sh ]]; then
  echo 'Please run this script from the top-level ExaGO source directory.'
  exit 1
fi

MAKEARGS=${MAKEARGS:-"-j"}
CTESTARGS=${CTESTARGS:-"--output-on-failure"}
EXTRA_CMAKE_ARGS=${EXTRA_CMAKE_ARGS:-""}
export OMPI_MCA_btl="^vader,tcp,openib,uct"
export BUILD=${BUILD:-1}
export TEST=${TEST:-1}
export CHECK_CMAKE=${CHECK_CMAKE:-0}
export SRCDIR=${SRCDIR:-$PWD}
export BUILDDIR=${BUILDDIR:-$PWD/build}
export INSTALLDIR=${INSTALLDIR:-$PWD/install}
export BUILD_MATRIX=${BUILD_MATRIX:-0}
export JOB=gcc-cuda
export VALID_JOBS=(gcc-cuda clang-hip clang-omp cmake-lint cmake-lint-apply)
export BUILD_VERBOSE=0

echo "Paths:"
echo "Source dir: $SRCDIR"
echo "Build dir: $BUILDDIR"
echo "Install dir: $INSTALLDIR"
echo "Path to buildsystem script: $SRCDIR/buildsystem/build.sh"
cd $SRCDIR

usage() {
  echo "Usage: ./buildsystem/build.sh [options]

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
  hostname deception.pnl.gov will be matched and 
  ./scripts/buildsystem/deceptionVariables.sh will be sourced. If you would like 
  to add a cluster, create a script
  ./scripts/buildsystem/<my cluster>Variables.sh and specify the relevant
  environment variables. If the hostname is not correctly finding your cluster,
  you may specify MY_CLUSTER environment variable before running this script
  and the script will respect the environment variable. For example, on ORNL
  Ascent cluster, the hostname does not find the cluster, so we must specify
  MY_CLUSTER when running:

    $ MY_CLUSTER=ascent ./buildsystem/build.sh --build-only

Spack:

  Each supported variables script in ./scripts/buildsystem activates a spack
  environment with all dependencies configured. If you have built dependencies
  for ExaGO in a spack environment, you may simply activate the environment
  and run the build script specifying that you don't want to source any
  variables scripts, eg:

    $ MY_CLUSTER=none ./buildsystem/build.sh

--------------------------------------------------------------------------------

Options:

  --job=<job name>  Run job indicated by job name. Available jobs are as
                    follows: ${VALID_JOBS[@]}.
                    Job --job=cmake-lint-apply should be ran before every push.

  --build-only      Only run the build stage of the script. This is useful for
                    local development.

  --test-only       Only run the test stage of the script. This should be ran
                    before every push to the repository or pull/merge request.
                    This run takes a significant amound of time. If you omit
                    the --*-only options and just run a particular job, tests
                    will also be ran.

  --verbose        Print all executed commands to the terminal. This is useful 
                   for debugging, but it will be disabled in CI by default to 
                   prevent hitting the job log limit.

--------------------------------------------------------------------------------

See ExaGO's latest developer guidelines for more information on developing
ExaGO: https://gitlab.pnnl.gov/exasgd/frameworks/exago/-/blob/develop/docs/DeveloperGuidelines.md

--------------------------------------------------------------------------------

"
}

while [[ $# -gt 0 ]]; do
  case $1 in
  --job=*)
    export JOB=$(echo $1 | cut -f2 -d'=')
    found=0
    for j in ${VALID_JOBS[@]}; do
      if [[ $JOB = $j ]]; then
        found=1
      fi
    done

    if [[ $found -eq 0 ]]; then
      echo
      echo "Job '$JOB' not available!"
      echo
      usage
      exit 1
    fi
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
  --help|-h)
    usage
    exit 0
    ;;
  --verbose|-v)
    export BUILD_VERBOSE=1
    shift
    ;;
  *)
    echo "Argument $1 not recognized."
    usage
    exit 1
    ;;
  esac
done

# Trap all signals and pass signal to the handler so we know what signal was
# sent in CI
for sig in `seq 1 31`; do
  trap "cleanup $sig \$? \$LINENO" $sig
done

if [ $BUILD_VERBOSE == 1 ]
then
  # verbose mode: print out all shell functions
  set -xv
else
  # don't print out all shell functions 
  set +x
fi

if [[ ! -v MY_CLUSTER ]]
then
  export MY_CLUSTER=`uname -n | sed -e 's/[0-9]//g' -e 's/\..*//'`
fi

# Correctly identify clusters based on hostname
# Valid cluster suggestions are spack based,
# however you can have any job with a valid:
# buildsystem/<job-name>/<platform>Variables.sh
case $MY_CLUSTER in
  newell*)
    export MY_CLUSTER=newell
    ;;
  incline*|dmi*)
    export MY_CLUSTER=incline
    ;;
  dl*|deception|*fat*)
    export MY_CLUSTER=deception
    ;;
  crusher*)
    export MY_CLUSTER=crusher
    ;; 
  ascent*)
    export MY_CLUSTER=ascent
    ;;
  summit*)
    export MY_CLUSTER=summit
    ;;
  *)
    echo "${MY_CLUSTER} did not match any directories in /buildsystem/spack/"
    echo "Try one of the following platforms: "
    echo $(ls -d ./buildsystem/spack/*/ | tr '\n' '\0' | xargs -0 -n 1 basename )
    exit 1 
    ;;
esac

ulimit -s unlimited || echo 'Could not set stack size to unlimited.'
ulimit -l unlimited || echo 'Could not set max locked memory to unlimited.'

varfile="$SRCDIR/buildsystem/$JOB/$(echo $MY_CLUSTER)Variables.sh"

if [[ -f "$varfile" ]]; then
  # source varfile without stderr or stout if it exists, error if failure
  set -xv
  source $varfile || { echo "Could not source $varfile"; exit 1; }
  # source $varfile 2>/dev/null || { echo "Could not source $varfile"; exit 1; }
  set +xv
fi

# module list

#  NOTE: The following is required when running from Gitlab CI via slurm job
base_path=`dirname $0`
if [ -z "$SLURM_SUBMIT_DIR" ]; then
  cd $base_path          || exit 1
fi

if [[ ! -f "$SRCDIR/buildsystem/$JOB/build.sh" ]]; then
  echo "Job $JOB not found!"
  usage
  exit 1
fi

source $SRCDIR/buildsystem/$JOB/build.sh
doBuild
EXIT_CODE=$?

if [ $BUILD_VERBOSE == 1 ]
then
  # turn off verbosity after build script 
  set +x
fi

exit $EXIT_CODE
