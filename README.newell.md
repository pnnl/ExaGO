# Installation and running ExaGO on Newell cluster

## System description

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

## Environment variables

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
libraries. To get needed modules
```console
$ module purge
$ module load gcc/7.4.0
$ module load openmpi/3.1.5
$ module load cmake
$ module load ipopt
```
Optimized PETSc library is installed on Newell in directory
```console
/qfs/projects/exasgd/newell/petsc
```
and debug version in
```console
/qfs/projects/exasgd/newell/petsc-dbg
```
You need to be member of `exasgd` group to access it. You can install your own
version of PETSc by downloading version >= 3.13 from [PETSc site](https://www.mcs.anl.gov/petsc/download/index.html)
and then configuring it like this:
```console
./configure                                    \
--with-cc=mpicc                                \
--with-cxx=mpicxx                              \
--with-fc=mpif90                               \
--with-debugging=0                             \
--with-cxx-dialect=C++11                       \
--download-fblaslapack                         \
--download-metis                               \
--download-mumps                               \
--download-parmetis                            \
--download-scalapack                           \
--download-suitesparse                         \
--download-superlu                             \
--download-superlu_dist                        \
--prefix=$INSTALL_DIR
```
where `$INSTALL_DIR` is a directory of your choosing. Once PETSc is configured,
just follow the instructions to compile, install, and test your PETSc build.

These instructions assume you loaded modules as described above.

## Install

ExaGO can be installed either on Newell easily with `cmake`.

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

The CMake installer should specify dynamic run paths. 


## Automatic setup upon login on Newell

To automate setting up development environment, you can add something like this
to your `.bashrc` shell startup script.

```bash
export MYSHORTHOST=`hostname -s | cut -c1-6`
if [ $MYSHORTHOST ==  "newell" ]; then
unset MODULE_VERSION
unset MODULE_VERSION_STACK
unset MODULESHOME
unset MODULEPATH
source /etc/profile.d/modules.sh
module purge
module load gcc
module load openmpi
export OMPI_MCA_pml="ucx"
export OMPI_MCA_btl="^vader,tcp,openib,uct"
export UCX_NET_DEVICES=mlx5_1:1,mlx5_3:1
fi
```

This will check if you are logged on Newell, load default GCC and OpenMPI
modules, and then set up environment variables needed by OpenMPI. You can edit
this script to meet your needs.

C-shell equivalent is here:
```shell
setenv MYSHORTHOST `hostname -s | cut -c1-6`
if ( $MYSHORTHOST =~  "newell" ) then
unsetenv MODULE_VERSION
unsetenv MODULE_VERSION_STACK
unsetenv MODULESHOME
unsetenv MODULEPATH
source /etc/profile.d/modules.csh
module purge
module load gcc
module load openmpi
setenv OMPI_MCA_pml "ucx"
setenv OMPI_MCA_btl "^vader,tcp,openib,uct"
setenv UCX_NET_DEVICES mlx5_1:1,mlx5_3:1
endif
```
