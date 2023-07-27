#!/bin/bash

#SBATCH -A EXASGD
#SBATCH -p incline
#SBATCH -n 3
#SBATCH --gres=gpu:1
#SBATCH --ntasks-per-gpu=1
#SBATCH --exclusive
#SBATCH -J exago_spack
#SBATCH -o spack_install.%J.output
#SBATCH -e spack_install.%J.output
#SBTACH -t 240

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

# Assuming that you already have a binary mirror configured (if you need it)
export MY_CLUSTER=incline
module purge
module load gcc/8.4.0
source buildsystem/clang-hip/inclineVariables.sh

# rocm-smi

# echo "incline variables loaded"
cd build
echo "========================================================"
id -a
echo "========================================================"
ctest -R UNIT_TESTS_OPFLOW_case9mod.m --verbose
# ctest -R FUNCTIONALITY_TEST_PFLOW_TESTSUITE_2_proc --verbose
# ctest -I 1
# ctest -I 15,16 # multinode fails
# ctest -I 5,29,33,34 failed tests

# MY_CLUSTER=incline buildsystem/build.sh --test-only --job=clang-hip

# printenv > srun_env.txt
whoami

EXIT_CODE=$?
# Required to trigger trap handler
exit $EXIT_CODE
