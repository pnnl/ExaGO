#!/bin/bash

#BSUB -P csc359
#BSUB -W 2:00
#BSUB -nnodes 1
#BSUB -J exasgd_spack_install
#BSUB -o spack_install.%J
#BSUB -e spack_install.%J

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

# Configure https proxy because spack is going to do some things with git
export all_proxy="socks://proxy.ccs.ornl.gov:3128"
export ftp_proxy="ftp://proxy.ccs.ornl.gov:3128"
export http_proxy="http://proxy.ccs.ornl.gov:3128"
export https_proxy="http://proxy.ccs.ornl.gov:3128"
export HTTP_PROXY="http://proxy.ccs.ornl.gov:3128"
export HTTPS_PROXY="http://proxy.ccs.ornl.gov:3128"
export proxy="proxy.ccs.ornl.gov:3128"
export no_proxy='localhost,127.0.0.0/8,*.ccs.ornl.gov,*.olcf.ornl.gov,*.ncrc.gov'


# Assuming that you already have a binary mirror configured
# TODO - copy over coinhsl tarball beforehand?
export MY_CLUSTER=summit
. buildsystem/spack/load_spack.sh && \
spack develop --no-clone --path=$(pwd) exago@develop && \
buildsystem/spack/configure_modules.sh 32

EXIT_CODE=$?
# Required to trigger trap handler
exit $EXIT_CODE
