# ExaGO (<b>Exa</b>scale <b>G</b>rid <b>O</b>ptimization)
ExaGO is a package for solving large-scale power grid optimization problems on parallel and distributed architectures, particularly targeted for exascale machines. It is wrritten in C/C++ using the [PETSc](https://www.mcs.anl.gov/petsc/) library. ExaGO includes the following power grid applications:

- [SC-ACOPF or SCOPFLOW](docs/web/scopflow.md) solves a security-constrained (contingency) constrained optimal power
- [ACOPF or OPFLOW](docs/web/opflow.md) solves an AC optimal power flow. OPFLOW is the basic building block for the SCOPFLOW application
- [ACPF or PFLOW](docs/web/pflow.md) solves an AC power flow

ExaGO applications can use the following solvers:

- [IPOPT](https://github.com/coin-or/Ipopt) is a popular optimization package for solving nonlinear optimization problems that uses an interior-point algorithm.
- [HiOp](https://github.com/LLNL/hiop) is a HPC solver for optimization. OPFLOW uses HiOp's mixed sparse-dense interior-point solver (NewtonMDS) that allows part of the computation to be run on GPU and part on CPU.
- [PIPS-NLP](https://github.com/Argonne-National-Laboratory/PIPS/tree/master/PIPS-NLP) is an HPC solver for structured optimization problems, such as those found in security-constrained optimization.
- [PETSc/TAO](https://www.mcs.anl.gov/petsc/) is a high-performance library providing numerical solvers for linear, nonlinear, time-dependent, and optimizatin problems. TAO includes a parallel primal-dual interior method (TAOPDIPM) that can be used for solving OPFLOW.


|  Solver | OPFLOW   | SCOPFLOW  | PFLOW  | Parallel  | GPU |
|:----------:|:---:|:---:|:---:|:---:|:---:|
| IPOPT      | Y   | Y   |     |     |   |
| HIOP       | Y   |     |     | Y   | Y |
| PIPS-NLP   |     | Y   |     | Y   |   |
| PETSc/TAO  | Y   |     | Y   | Y   | Y |

## Prerequisites
ExaGO depends on various third-party packages. PETSc is a core-dependency of ExaGO. At least one of the solver packages (IPOPT, HiOp, or PIPS-NLP) should be installed depending the application.
1. [PETSc ver. 3.13](docs/web/petsc_install.md)
1. [IPOPT](docs/web/ipopt_install.md) 
1. [HiOp](docs/web/hiop_install.md)
1. [PIPS-NLP](docs/web/pips_install.md)

In addition, we recommend installing `git` and `cmake` (ver. 3.10 or higher)

## Download
Download ExaGO from gitlab via
```
git clone https://gitlab.pnnl.gov/exasgd/frameworks/exago.git
```
or if you have SSH access to the repository
```
git clone ssh://git@gitlab.pnnl.gov:2222/exasgd/frameworks/exago.git
```

## Install

ExaGO can be installed either with `cmake` (preferred) or `make`.

### Building with `cmake`

Note: CMake instructions are work in progress and will continue to be updated

First create build directory outside the ExaGO source directory. For example
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
### Building with `make`
Another way to build ExaGO is to use `make`. This is purely for development purposes and will be replace with `cmake` once it is stable. For building with `make`, an environment variable for ExaGO directory location needs to be set.
```
export EXAGO_DIR=<location of ExaGO>
```
Environment variables for [PETSc](docs/web/petsc_install.md), [IPOPT](docs/web/ipopt_install.md), [HIOP](docs/web/hiop_install.md), and [PIPS](docs/web/pips_install.md). 

Individual application can then be compiled with the `make` command.
```
cd $EXAGO_DIR
make <appliction_name>
```
Appication name is either `ExaGO`, `OPFLOW`, or `PFLOW`. Use flags `WITH_IPOPT=1`,`WITH_HIOP=1`, or `WITH_PIPS=1` to enable IPOPT, HiOp, and PIPS, respectively. For example,
```
make OPFLOW WITH_HIOP=1
```
will build OPFLOW application enabling HiOp package.

#### Additional installation tips (Linking third-party libraries)
Depending on how you installed PETSc, PIPS and IPOPT, you may need to update the linker environment variables to include the PETSc, IPOPT, PIPS-NLP, and ExaGO library paths
```
export LD_LIBRARY_PATH=$PETSC_DIR/$PETSC_ARCH/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$EXAGO_DIR:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$IPOPT_BUILD_DIR/lib:$LD_LIBRARY_PATH
LD_LIBRARY_PATH=$PIPS_DIR/build/PIPS-NLP:$LD_LIBRARY_PATH
```

## Usage

- [SCOPFLOW](docs/web/scopflow.md)
- [OPFLOW](docs/web/opflow.md)
- [PFLOW](docs/web/pflow.md)

## Authors

## Acknowledgement
This package is developed as a part of [ExaSGD](https://www.exascaleproject.org/wp-content/uploads/2019/10/ExaSGD.pdf) project under the [Exascale computing project](https://www.exascaleproject.org/).

## License


