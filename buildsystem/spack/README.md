# Spack Build Scripts

Run all scripts from root git directory.

If you are installing with CoinHSL, ensure you have applied for a personal/academic
license. Then, copy the CoinHSL tarball into the ExaGO root directory.

You will have to re-name this tarball to match the name that spack expects, and
you can find this name by running `spack edit coinhsl` once spack is loaded.

These scripts work using the spack submodule in `tpl/spack`. This maintains
a reference commit of spack so that one can re-produce builds, and share
installations on shared file systems. Make sure to recursively clone
submodules before running using `git submodule update --init --recursive`.

We would just use a tagged version of spack, however ExaGO undergoes rapid
development alongside HiOp, and so this is the approach that gives the most
flexibility.

Whenever an update to the ExaGO spack submodule reference is made, it may render
any existing installations unusable due to changes in the hashing algorithm.

## Using the .pre pipelines in GitLab PRs

Whenever a new PR is created, you can manually run a set of GitLab pipelines that will
automatically generate updated spack modules and push changes to
the branch that triggered the pipeline to use the newly generated modules. This commit 
will trigger new pipelines on each platform, so you may want to manually cancel CI on
platforms that were already tested in the pipeline that created the new modules.

- TODO: add commit messages to each module update to only run pipeline on updated platform

## General Workflow

Run scripts in the following order to install using spack for a supported platform.

**Make sure to run `export MY_CLUSTER=platform` before running scripts.**

1. `load_spack.sh`

This should be run first, and will create a spack environment in the git root
ExaGO directory under `spack-env-platform`.

This pulls spack config from `buildsystem/spack/platform/spack.yaml` into your
environment, and installs spack based on environment variables loaded in
`buildsystem/spack/platform/env.sh`.


2. `spack concretize -f`

This is a regular spack command, and should be run to verify what configuration
is being built, before committing to an installation.


4. `configure_modules.sh`

This will build the config of ExaGO with and without optimizations, and then
populate the related shell modules with the relevant module path, and necessary
module commands.


## Installing on a compute node/without internet access

1. `binary_mirror.sh` 

This will configure a binary mirror based on the current environment, and is
required before installing on a resource without internet connection (i.e. a 
compute node in a secure enclave).

The environment variable `SOURCE_CACHE` will specify the specific folder in
`SPACK_INSTALL` that the binary mirror will be configured for. This should
be set in the `env.sh` script for your platform.

2. `job script` 

There should be a job script that you can submit to the job scheduler that will
automatically install on a compute node using spack, given you have configured
a binary mirror where necessary.

This will vary depending on the platforms job scheduler, and there will be a
custom config for each platform.

3. `install_compute.sh` 

Since we often need to install using spack on a compute node, this script
automates that process for each platform. Make sure to verify your environment
is correct before running this, and that a binary mirror is configured is running
on an enclave without internet access.

Note: This calls the job script for a given platform, which will build
with `spack develop`. This allows for a module set generation with binaries
built from the current branch in the git clone of ExaGO.

## Configuring for a new platform

In order to add configuration for a new platform, you need to add the following:

- `env.sh`

    - All configurations for a platform that can be controlled in environment variables 
    is in the `buildsystem/spack/platform/env.sh` script.

- `spack.yaml`

    - A fully descriptive environment file that stores compiler information, along 
    with specification of the ExaGO builds and dependencies for a given platform.

- `install.sh`

    - A simple script that acts as an intermediary for submitting jobs. You can
     just submit the job script yourself - this is a convenience when dealing
     with many different queueing systems.

- `job_script`

    - The full `sbatch` (or `bsub`) script used to install modules. This often is
     configured to have specific error handling if tied to a CI job, along with 
     using `spack develop` to install based on the current branch.

Once you have all of this for a new platform, you will get various module files
generated in `buildsystem/spack/platform/modules`. This includes 4 module
configurations currently: exago dependencies, exago standalone, optimized exago
dependencies, and optimized exago standalone. This is dependent on the configuration
in `spack.yaml` and you should make sure to have a consistent spack command
in `configure_modules.sh` that aligns with that is specified in the environment.

These modules are used in `buildsystem/<build_type>/` in various scripts for each
platform, and you should ensure these are compatible with the module configuration
used for each platform. Make sure to create variable files that are consistent
with examples for `gcc-cuda/newell` and `clang-hip/crusher`.

If you have an update in the ExaGO/HiOp spack package, you might need to update
the relevant `spack.yaml` and configure module scripts for each platform.

