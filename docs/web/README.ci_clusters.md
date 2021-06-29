<!--- vi: set filetype=markdown: -->

# Introduction

This document outlines building and installing ExaGO on three clusters:

- Newell (Power9 PNNL system)
- Marianas (Intel PNNL system)
- Ascent (Power9 ORNL system)

This document goes into great detail about building on Newell, but most of the
instructions apply to all clusters.

# Quick Start

If you are on one of the aforementioned systems, the following commands should
yield an ExaGO installation with most or all options enabled, including GPU
computation.

```console
$ # Set this variable to one of newell, marianas, or ascent
$ export MY_CLUSTER=newell

$ git clone https://gitlab.pnnl.gov/exasgd/frameworks/exago.git
$ cd exago
$ mkdir build install

$ # Load all the modules needed to build/run ExaGO
$ source ./buildsystem/gcc-cuda/${MY_CLUSTER}Variables.sh
$ cd build

$ # Use the initial CMake cache we use for CI
$ cmake -C ../buildsystem/gcc-cuda/cache.cmake ..

$ make -j 12 install

$ # The tests may take a while to run
$ make test
```

# Additional Information

## Install and run ExaGO on Newell cluster

### System description

Newell is a cluster of five IBM AC922 Power9 servers four NVIDIA Volta GPUs per
system. This system gives researchers unfettered access to
architectures found in the Summit supercomputer at Oak Ridge that is currently
ranked No. 1 on the [Top500 list](https://www.top500.org/).

Specs:

- 2 Power9 CPUs with a total of 128 logical cores per node
- 4 NVIDIA V100 GPUs with NVLINK (16GB per GPU)
- 1TB of system memory per node
- EDR Infiniband internal network
- 10Gb/s connections to PNNL network
- Nodes have 2 physical Infiniband ports each of which has 2 virtual ports

More information is available at Research Computing
[confluence page](https://confluence.pnnl.gov/confluence/display/RC/Newell).

### Quick start using CI Configuration

More details on setting up your environment to run on the Newell cluster are
available below, but the current spack packages can be used to build ExaGO and run
the ExaGO test suite by typing just a few lines.

The build scripts currently only work in the bash shell. If you are not using
bash, type `bash` immediately after logging into newell. Then `cd` into the
root source directory of ExaGO.

Then source the shell script we use to load the needed modules and set
the needed environment variables to build on Newell
```
source buildsystem/gcc-cuda/newellVariables.sh
```
This will set all environment variables needed by ExaGO.

Next, request an interactive session on one of the compute nodes:
```
srun -A exasgd -t 20:00 --gres=gpu:1 -p newell -n 2 --pty bash
```

The <code>-A exasgd</code> indicates you are using the ExaSGD allocation, <code>-t
20:00 </code> specifies that you are requesting 20 minutes for your interactive
session, <code>--gres=gpu:1</code> configures the GPUs, <code>-p newell</code> means
that you are using the newell partition on SLURM, <code>-n 2</code> specifies the
number of MPI tasks for this session that you are requesting and
<code>--pty bash</code> means your interactive session will be using the bash shell.
One or two tests are using two MPI task, so you sould select at least two tasks
in the interactive shell.  Using the newell partition guarantees that all your
environment variables and binary files match up with the hardware you are running on.

To build ExaGO, just type
```
./buildsystem/build.sh --job=gcc-cuda --build-only
```
in the top level directory while in the interactive session. This will create a
<code>build</code> directory underneath the ExaGO directory. This directory contains
all program executables, makefiles and test directories. In addition to
building ExaGO, the script will also run the test suite. If you wish to rerun
the test suite after running the <code>build.sh</code> script you can do so by
ommitting the `--build-only` argument to the build script, or by going
into the <code>build</code> directory and typing either
```
ctest
```
or
```
make test
```

### Environment variables

Currently, runing OpenMPI on Newell requires some nonstandard flags to be
passed to the `mpirun` script. The MPI run command
```
mpirun -n 4
```
may need to be replaced with something like this:
```
mpirun -mca pml ucx --mca btl ^vader,tcp,openib,uct  -x UCX_NET_DEVICES=mlx5_1:1,mlx5_3:1 -n 4
```
Without setting these flags, OpenMPI might report warnings or even fail when
running multi-node jobs.

You can save yourself some typing by setting environment variables like this
```
export OMPI_MCA_pml="ucx"
export OMPI_MCA_btl="^vader,tcp,openib,uct"
export UCX_NET_DEVICES=mlx5_1:1,mlx5_3:1
```

For more information take a look at OpenMPI [FAQs](https://www.open-mpi.org/faq/?category=tuning#mca-def).

## Modules

Currently tested (and recommended) tool chain for building ExaGO consists of
GCC 7.4 and OpenMPI 3.1.5. You will also need recent version of CMake. ExaGO
depends on PETSc >= 3.13, and optionally on Ipopt and HiOp optimization
libraries. To get needed modules, refer to the CI scripts under `buildsystem/gcc-cuda/`.

## Install

ExaGO can be installed on Newell easily with `cmake`.

First create build directory outside the ExaGO source directory. For example
```
$ mkdir build
$ ls
build  exago
$
```
Then from build directory configure ExaGO using `cmake`:
```
$ cd build
$ cmake ../exago
$ make install
```

If you would like to configure ExaGO to use all options used in continuous
integration, invoke CMake like so:

```console
$ cd build
$ cmake ../exago -C ../exago/buildsystem/gcc-cuda/cache.cmake
$ make install
```

The ExaGO library and its applications are installed in the default
installation directory. To change installation directory run CMake with flag
```
$ cmake ../exago -DCMAKE_INSTALL_PREFIX=<your_exago_install_dir>
```
ExaGO assumes PETSc is built with MPI support. If it is not, it is recommended
you configure ExaGO not to use MPI: 
```
$ cmake -DEXAGO_ENABLE_MPI=Off ../exago
```
In case PETSc dependency is not automatically found, you can specify it using
`ccmake` interactive shell or add command line option like this:
```
$ cmake ../exago -DPETSC_DIR=<petsc_install_dir> -DPETSC_ARCH=<petsc_arch>
```

To use IPOPT with ExaGO, set:
```
cmake ../exago -DEXAGO_ENABLE_IPOPT=ON
```
ExaGO will find Ipopt module you loaded on Newell. If you want to use your
own Ipopt build, you will most likely need to specify its location like this:
```
cmake ../exago -DEXAGO_ENABLE_IPOPT=ON -DIPOPT_DIR=<ipopt_install_dir>
```
Similar to IPOPT, the corresponding flags for HiOp are `EXAGO_ENABLE_HIOP`
and `HIOP_DIR`.
```
cmake ../exago -DEXAGO_ENABLE_HIOP=ON -DHIOP_DIR=<hiop_install_dir>
```

The CMake installer should specify dynamic run paths. To run ExaGO, you may need
to load all the modules you used to build it.

## Notes on Ascent

Ascent uses the LSF scheduler for job submission, so if you would like to run
the tests, you'll have to pass some additional options.

```console
$ export MY_CLUSTER=ascent

$ git clone https://gitlab.pnnl.gov/exasgd/frameworks/exago.git
$ mkdir build install
$ cd exago

$ # Load all the modules needed to build/run ExaGO
$ source ./buildsystem/gcc-cuda/${MY_CLUSTER}Variables.sh
$ cd ../build

# The EXAGO_TEST_WITH_BSUB option ensures that `make test` will use `jsrun`
# to launch each test on a compute node.
$ cmake \
  -C ../exago/buildsystem/gcc-cuda/cache.cmake \
  -DEXAGO_TEST_WITH_BSUB=ON \
  -DCMAKE_INSTALL_PREFIX=$PWD/../install \
  ../exago

$ make -j 12 install

$ # Request an allocation for 15 minutes
$ bsub -P csc359 -W 15 -nnodes 1 -Is /bin/bash

$ # The tests may take a while to run
$ make test
```

## Notes on Marianas

A workflow on Marianas should look almost exactly like on Newell, except you
should request an allocation in the `dl` or `dl_shared` partition, and
source `buildsystem/gcc-cuda/marianasVariables.sh`.
