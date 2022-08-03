## ExaGO Buildsystem 

This directory contains tools used to manage the deployment of the ExaGO codebase.

This documentation is under development. Some brief notes will be added about what platform 
is supported by each configuration, or a short description of what the purpose is.

### clang-omp

Platforms:

- Newell

Description:

newell clang build of exago@1.0.0 + hiop@0.4.0

### cmake

Description:

Contains any CMake configuration for configuring external libraries.

### container

Platforms:

Used in spack-ci stage of testing (currently not functioning).

Description:

https://gitlab.pnnl.gov/exasgd/frameworks/exago/-/commit/47ea09e648dfa81ca8a70cc2e0b4a2628eaf56f2 - here is the merge request into master that created this initially. Remnant of previously functional spack-ci.

### gcc-cuda

Platforms:
- Ascent
- Deception
- Marianas
- Newell
- Summit

Description:

gcc and cuda enabled build for target x86_64 and power9 platforms.

### misc

Description:

Old miscellaneous scripts

### pnnl

Description:

Internal perl scripts that we use to find idle nodes within our clusters using slurm. Are a component of CI pipelines on Marianas and Newell.

### tools - ExaGO Code Quality Tools

This directory contains tools used to manage the code quality of the ExaGO codebase. See README.md inside for more.

### build.sh

A self-documented script that can enable building on most target platforms:

```bash
$ ./buildsystem/build.sh --help
Paths:
Source dir: /ccs/home/rcruther/exago-git
Build dir: /ccs/home/rcruther/exago-git/build
Install dir: /ccs/home/rcruther/exago-git/install
Path to buildsystem script: /ccs/home/rcruther/exago-git/buildsystem/build.sh
Usage: ./buildsystem/build.sh [options]

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
                    follows: gcc-cuda clang-hip clang-omp cmake-lint cmake-lint-apply.
                    Job --job=cmake-lint-apply should be ran before every push.

  --build-only      Only run the build stage of the script. This is useful for
                    local development.

  --test-only       Only run the test stage of the script. This should be ran
                    before every push to the repository or pull/merge request.
                    This run takes a significant amound of time. If you omit
                    the --*-only options and just run a particular job, tests
                    will also be ran.

--------------------------------------------------------------------------------

See ExaGO's latest developer guidelines for more information on developing
ExaGO: https://gitlab.pnnl.gov/exasgd/frameworks/exago/-/blob/develop/docs/DeveloperGuidelines.md

--------------------------------------------------------------------------------


BUILD_STATUS:0
```
