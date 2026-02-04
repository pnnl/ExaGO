## ExaGO Buildsystem 

This directory contains tools used to manage the deployment of the ExaGO codebase.

This documentation is under development. Some brief notes will be added about what platform 
is supported by each configuration, or a short description of what the purpose is.

Each folder which builds a configuration of ExaGO should have a following:
- `build.sh` : custom build script for configuration. Usually template
- `platformVariables.sh` : a list of spack modules and other config to set up a CMake build environment.
- `cache.cmake` : universal CMake changes for the type of build in use
- `platform` :
  - `base.sh` : a script to load platform-specific modules
  - `platformExago.sh` : a script to load that build of exago for use directly on command line

### clang-hip

Platforms:

- frontier

Description:

frontier clang build of exago@frontier-dev + hiop@develop

### cmake

Description:

Contains any CMake configuration for configuring external libraries.

### container

Platforms:

Used in spack-ci stage of testing (currently not functioning).

### deprecated

Description:

Contains deprecated build configurations.

#### clang-omp

Platforms:

- Newell

Description:

Newell clang build of exago@1.0.0 + hiop@0.4.0

#### gcc-cuda

Platforms:
- Deception
- Marianas
- Newell

Description:

gcc and cuda enabled build for target x86_64 and power9 platforms.

### misc

Description:

Old miscellaneous scripts

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
  variable MY_CLUSTER will be set. If you would like to add a cluster, create a script
  ./buildsystem/<job name>/<my cluster>Variables.sh and specify the relevant
  environment variables. If the hostname is not correctly finding your cluster,
  you may specify MY_CLUSTER environment variable before running this script
  and the script will respect the environment variable. For example, on ORNL
  Frontier cluster, the hostname does not find the cluster, so we must specify
  MY_CLUSTER when running:

    $ MY_CLUSTER=frontier ./buildsystem/build.sh --build-only

Spack:

  Each supported variables script in ./buildsystem/<job name> activates a spack
  environment with all dependencies configured. If you have built dependencies
  for ExaGO in a spack environment, you may simply activate the environment
  and run the build script specifying that you don't want to source any
  variables scripts, eg:

    $ MY_CLUSTER=none ./buildsystem/build.sh

--------------------------------------------------------------------------------

Options:

  --job=<job name>  Run job indicated by job name. Available jobs are as
                    follows: clang-hip.
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


BUILD_STATUS:0
```
