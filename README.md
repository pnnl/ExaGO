# <b>Exa</b>scale <b>G</b>rid <b>O</b>ptimization toolkit (ExaGO<sup>TM</sup>)
ExaGO<sup>TM</sup> is a package for solving large-scale power grid optimization problems on parallel and distributed architectures, particularly targeted for exascale machines with heteregenous architectures (GPU). It is written in C/C++ using the [PETSc](https://www.mcs.anl.gov/petsc/) library. ExaGO<sup>TM</sup> includes the following power grid applications:

- [S-ACOPF or SOPFLOW](docs/web/sopflow.md) solves a stochastic security-constrained multi-period optimal power flow
- [SC-ACOPF or SCOPFLOW](docs/web/scopflow.md) solves a multi-period security-constrained (contingency) constrained optimal power
- [TC-ACOPF or TCOPFLOW](docs/web/tcopflow.md) solves a multi-period optimal power flow
- [ACOPF or OPFLOW](docs/web/opflow.md) solves an AC optimal power flow either on CPU and GPU
- [ACPF or PFLOW](docs/web/pflow.md) solves an AC power flow



ExaGO<sup>TM</sup> applications can use the following solvers:

- [IPOPT](https://github.com/coin-or/Ipopt) is a popular optimization package for solving nonlinear optimization problems that uses an interior-point algorithm.
- [HiOp](https://github.com/LLNL/hiop) is a HPC solver for optimization. OPFLOW uses HiOp's mixed sparse-dense interior-point solver (NewtonMDS) that allows part of the computation to be run on GPU and part on CPU.
- [PETSc/TAO](https://www.mcs.anl.gov/petsc/) is a high-performance library providing numerical solvers for linear, nonlinear, time-dependent, and optimizatin problems. The optimization package is PETSc, TAO, includes a parallel primal-dual interior method (TAOPDIPM) that can be used for solving OPFLOW.


|  Solver    | OPFLOW | SCOPFLOW  | PFLOW | TCOPLOW | SOPFLOW |
|:----------:|:------:|:---------:|:-----:|:-------:|:-------:|
| IPOPT      | Y      | Y         |       | Y       | Y       | 
| HIOP       | Y      |           |       |         |         |
| PETSc/TAO  | Y      |           | Y     |         |         |


## Prerequisites
ExaGO depends on the [PETSc](docs/web/petsc_install.md) library and is compatible with version 3.13.

In addtion, at least one of the solver packages should be installed depending the application.
1. [IPOPT](docs/web/ipopt_install.md) (version 3.12 and above)
1. [HiOP](docs/web/hiop_install.md) (develop branch)

ExaGO can run the OPFLOW application on the GPU using HiOP solver library and, in addition, needs UMPIRE and RAJA libraries that provide portability layer and memory management for running calculations on the GPU.
1. [RAJA](https://github.com/LLNL/RAJA)
1. [UMPIRE](https://github.com/LLNL/Umpire) 

## Download
Download ExaGO<sup>TM</sup> from gitlab via
```
git clone https://gitlab.pnnl.gov/exasgd/frameworks/exago.git
```
or if you have SSH access to the repository
```
git clone ssh://git@gitlab.pnnl.gov:2222/exasgd/frameworks/exago.git
```

## Building and Installing with CMake
ExaGO uses a CMake build system for building, installing, and testing. We recommend using `cmake` ver. 3.10 or higher. To build
ExaGO with CMake, first create build directory outside the ExaGO source directory. For example
```shell
$ mkdir build
$ ls
build  exago
$
```
Then from build directory configure ExaGO using `cmake`:
```shell
$ cd build
$ cmake ../exago
$ make install
```
The ExaGO library and its applications are installed in the default installation
directory. To change installation directory run CMake with flag
```
$ cmake ../exago -DCMAKE_INSTALL_PREFIX=<your_exago_install_dir>
```
ExaGO assumes PETSc is built with MPI support. If it is not, it is recommended
you configure ExaGO not to use MPI: 
```
$ cmake -DEXAGO_ENABLE_MPI=Off ../exago
```

To use ExaGO without MPI, you must also build PETSc without MPI. See [PETSc installation](docs/web/petsc_install.md) for instructions on how to build PETSc without MPI.


In case PETSc dependency is not automatically found, you can specify it using
`ccmake` interactive shell or add command line option like this:
```
$ cmake ../exago -DPETSC_DIR=<petsc_install_dir> -DPETSC_ARCH=<petsc_arch>
```

To use IPOPT with ExaGO, set the `EXAGO_ENABLE_IPOPT`. When this flag is set, ExaGO will try to find IPOPT in some default locations and will error if IPOPT is not found. In this case, the IPOPT installation directory should be set with `IPOPT_DIR'.
```
cmake ../exago -DEXAGO_ENABLE_IPOPT=ON -DIPOPT_DIR=<ipopt_install_dir>
```

Similar to IPOPT, the corresponding flags for HiOp are 'EXAGO_ENABLE_HIOP' and 'HIOP_DIR'.
```
cmake ../exago -DEXAGO_ENABLE_HIOP=ON -DHIOP_DIR=<hiop_install_dir>
```

Below is an example build with all (optional) dependencies installed
```
cmake ../
-DCMAKE_INSTALL_PREFIX=$installdir/ \
  -DCMAKE_BUILD_TYPE=Debug \
  -DEXAGO_ENABLE_GPU=ON \
  -DEXAGO_ENABLE_HIOP=ON \
  -DEXAGO_ENABLE_IPOPT=ON \
  -DEXAGO_ENABLE_MPI=ON \
  -DEXAGO_ENABLE_PETSC=ON \
  -DEXAGO_RUN_TESTS=ON \
  -DEXAGO_ENABLE_RAJA=ON \
  -DRAJA_DIR=$raja_dir \
  -Dumpire_DIR=$umpire_dir \
  -DHIOP_DIR=$hiop_dir \
  -DIPOPT_DIR=$ipopt_dir \
  -DMAGMA_DIR=$magma_dir \
  -DPETSC_DIR=$petsc_dir
```

## Usage
Instructions for executing the different ExaGO applications is given below.
- [OPFLOW](docs/web/opflow.md)
- [SOPFLOW](docs/web/sopflow.md)
- [SCOPFLOW](docs/web/scopflow.md)
- [PFLOW](docs/web/pflow.md)

## Authors
- Shrirang Abhyankar
- Slaven Peles
- Asher Mancinelli
- Robert Rutherford
- Bruce Palmer

## Acknowledgement
This package is developed as a part of [ExaSGD](https://www.exascaleproject.org/wp-content/uploads/2019/10/ExaSGD.pdf) project under the [Exascale computing project](https://www.exascaleproject.org/).

## Copyright
Copyright &copy; 2020, Battelle Memorial Institute.

ExaGO is a free software distributed under a BSD 2-clause license. You may reuse, modify, and redistribute the software. See the [license](LICENSE) file for details.


## Disclaimer
This material was prepared as an account of work sponsored by an agency of the United States Government.  Neither the United States Government nor the United States Department of Energy, nor Battelle, nor any of their employees, nor any jurisdiction or organization that has cooperated in the development of these materials, makes any warranty, express or implied, or assumes any legal liability or responsibility for the accuracy, completeness, or usefulness or any information, apparatus, product, software, or process disclosed, or represents that its use would not infringe privately owned rights.
Reference herein to any specific commercial product, process, or service by trade name, trademark, manufacturer, or otherwise does not necessarily constitute or imply its endorsement, recommendation, or favoring by the United States Government or any agency thereof, or Battelle Memorial Institute. The views and opinions of authors expressed herein do not necessarily state or reflect those of the United States Government or any agency thereof.
