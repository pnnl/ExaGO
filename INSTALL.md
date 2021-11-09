
# Installing ExaGO

## Table of Contents

1. [Installing With Spack](#installing-with-spack)
1. [Installing From Source](#installing-from-source)
    1. [Dependencies](#dependencies)
        1. [General dependencies](#general-dependencies)
        1. [Third-Party Solvers](#third-party-solvers)
        1. [GPU-related dependencies](#gpu-related-dependencies)
    1. [Additional CMake Options](#additional-cmake-options)
    1. [CMake Workflow](#cmake-workflow)
    1. [CMake Configuration Options](#cmake-configuration-options)
1. [Building on CI Platforms](#building-on-ci-platforms)
    1. [Additional Information](#additional-information)

## Installing With Spack

If you would like to use the Spack package manager to install ExaGO, please refer to [the ExaGO Spack documentation](./docs/InstallingWithSpack.md).
The rest of this document assumes you want to install with CMake and you already have the dependencies described in the [Dependencies section linked here](#dependencies).

## Installing From Source

### Acquiring the Source Code

```console
git clone https://gitlab.pnnl.gov/exasgd/frameworks/exago.git

# Or if you have SSH keys set up:
git clone ssh://git@gitlab.pnnl.gov:2222/exasgd/frameworks/exago.git
```

### Dependencies

Dependencies are broken down into three categories:

1. [General dependencies](#general-dependencies)
1. [Third-Party Solvers](#third-party-solvers)
1. [GPU-related dependencies](#gpu-related-dependencies)

#### General Dependencies

| Dependency | Dependency Version | CMake Variable (if applicable) | Notes |
|---|---|---|---|
| PETSC | >=3.13, <3.15 | `EXAGO_ENABLE_PETSC`. `PETSC_DIR` sets the PETSC installation prefix. | Required. PETSC needs to be built with MPI if building ExaGO with MPI. [See the linked document for more information on building PETSC for ExaGO.](docs/web/petsc_install.md) |
| CMake | >=3.18 | N/A |  |
| blas |  | `BLAS_LIBRARIES` sets the blas libraries. | ExaGO uses the builtin [`FindBLAS` cmake module](https://cmake.org/cmake/help/latest/module/FindBLAS.html). Please refer to the linked documentation. |
| MPI | >=3 | `EXAGO_ENABLE_MPI`. `MPI_HOME` sets the MPI installation prefix. | Optional. ExaGO uses the builtin [`FindMPI` cmake module](https://cmake.org/cmake/help/latest/module/FindMPI.html). |
| Python | >=3.5 | `EXAGO_ENABLE_PYTHON`, `Python3_DIR` | ExaGO uses builtin [`FindPython` cmake module](https://cmake.org/cmake/help/latest/module/FindPython.html). Please refer to linked documentation. |

#### Third-Party Solvers

At least one of the solver packages should be installed depending the application.
Additional documentation on installing these dependencies are linked below.

| Dependency | Dependency Version | CMake Variable (if applicable) | Notes |
|---|---|---|---|
| HiOp | v0.5 | `EXAGO_ENABLE_HIOP`, `HIOP_DIR` | Optional, needs to be built with CUDA and MPI if using those options. [See HiOp repo linked here for more information.](https://github.com/LLNL/hiop) |
| Ipopt | >=3.12 | `EXAGO_ENABLE_IPOPT`, `IPOPT_DIR` | Needs to be built with COINHSL if using HiOp sparse linear algebra. [See the linked document for more information on building Ipopt for ExaGO](./docs/web/ipopt_install.md), [and see this link for PETSC documentation.](https://petsc.org/release/)|

#### GPU-Related Dependencies

ExaGO can run the OPFLOW application on a GPU using the HiOP solver library and, in addition, needs UMPIRE and RAJA libraries that provide a portability layer over execution and memory management for running calculations on the GPU.
If you are using the HiOp solver, HiOp must be built with the same RAJA and Umpire options, and must use the same GPU backend.

| Dependency | Dependency Version | CMake Variable | Notes |
|---|---|---|---|
| CUDA | >=10.2 | `EXAGO_ENABLE_GPU` and `EXAGO_ENABLE_CUDA` | Optional. `EXAGO_ENABLE_GPU` and `EXAGO_ENABLE_CUDA` must both be enabled because other GPU platforms are under development. `CUDAToolkit_ROOT`  and `CMAKE_CUDA_COMPILER` can be used to select a specific CUDA installation if you have multiple. |
| RAJA | >=0.14 | `EXAGO_ENABLE_RAJA`, `RAJA_DIR` | Optional, needs to be built with CUDA if using `EXAGO_ENABLE_CUDA`. [Refer to RAJA documentation linked here.](https://github.com/LLNL/RAJA) |
| Umpire | >=6.0 | `EXAGO_ENABLE_RAJA`, `umpire_DIR` | Optional, needs to be built with CUDA if using `EXAGO_ENABLE_CUDA`. [Refer to Umpire documentation linked here.](https://github.com/LLNL/Umpire) |
| MAGMA | >=2.6.1 | `EXAGO_ENABLE_HIOP`, `MAGMA_DIR` | Only required when using HiOp solve (`EXAGO_ENABLE_HIOP`) with `EXAGO_ENABLE_GPU`. Must be same version your HiOp was built with. |

### Additional CMake options

| CMake Variable | Notes |
|---|---|
| `CMAKE_INSTALL_PREFIX` | Standard CMake variable. Sets ExaGO installation prefix. |
| `CMAKE_BUILD_TYPE` | Standard CMake variable. Sets ExaGO build type, eg `Release`, `Debug`, or `RelWithDebInfo`. |
| `EXAGO_BUILD_SHARED` | Build ExaGO shared libraries. At least one of `EXAGO_BUILD_SHARED` and `EXAGO_BUILD_STATIC` must be enabled. |
| `EXAGO_BUILD_STATIC` | Build ExaGO static libraries. |
| `EXAGO_TEST_WITH_BSUB` | Use BSUB/JSRUN commands in CTest. Use if you would like to run ExaGO test suite in a BSUB allocation. |
| `EXAGO_EXTRA_MPI_FLAGS` | Pass extra flags to `mpirun`/`mpiexec`/`jsrun`/`srun` in CTest.  |

### CMake Workflow

```console
git clone ssh://git@gitlab.pnnl.gov:2222/exasgd/frameworks/exago.git
cd exago
mkdir build
cd build

# Pass any configuration options here
cmake .. -DEXAGO_BUILD_SHARED=ON

# You may also cusomize your build here
ccmake . 

make -j
make install
make test # if you like
```

### CMake Configuration Options

ExaGO uses a CMake build system for building, installing, and testing.
To build ExaGO with CMake, first create build directory outside the ExaGO source directory. For example
```console
mkdir build
```
Then from build directory configure ExaGO using `cmake`:
```shell
cd build
cmake ../exago
make install
```

The ExaGO library and its applications are installed in the default installation
directory. To change installation directory run CMake with flag
```console
$ cmake ../exago -DCMAKE_INSTALL_PREFIX=<your_exago_install_dir>
```

ExaGO assumes PETSc is built with MPI support. If it is not, it is recommended
you configure ExaGO not to use MPI: 
```console
$ cmake -DEXAGO_ENABLE_MPI=Off ../exago
```

To use ExaGO without MPI, you must also build PETSc without MPI. See [PETSc installation](docs/web/petsc_install.md) for instructions on how to build PETSc without MPI.

In case PETSc dependency is not automatically found, you can specify it using
`ccmake` interactive shell or add command line option like this:
```console
$ cmake ../exago -DPETSC_DIR=<petsc_install_dir> -DPETSC_ARCH=<petsc_arch>
```

To use IPOPT with ExaGO, set the option `EXAGO_ENABLE_IPOPT`. When this flag is set, ExaGO will try to find IPOPT in some default locations and will error if IPOPT is not found. In this case, the IPOPT installation directory should be set with `IPOPT_DIR'.
```console
cmake ../exago -DEXAGO_ENABLE_IPOPT=ON -DIPOPT_DIR=<ipopt_install_dir>
```

Similar to IPOPT, the corresponding flags for HiOp are 'EXAGO_ENABLE_HIOP' and 'HIOP_DIR'.
```console
cmake ../exago -DEXAGO_ENABLE_HIOP=ON -DHIOP_DIR=<hiop_install_dir>
```

Below is an example build with all (optional) dependencies installed
```console
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

## Building on CI Platforms

***NOTE: For developers only***

If you are building on one of our continuous integration (CI) platforms, you may
want to use the CI environment for development. To do this, source the variables
script like so:

```console
$ # On newell for example
$ source ./buildsystem/gcc-cuda/newellVariables.sh
$ mkdir build && cd build && cmake ..
```

If you would like to use the *exact* configuration used in CI, you may use the
cmake cache script like so:

```console
$ # On newell for example
$ source ./buildsystem/gcc-cuda/newellVariables.sh
$ mkdir build && cd build
$ cmake .. -C ../buildsystem/gcc-cuda/cache.cmake
```

However, users are encouraged to configure their build on their own. Please use
this CMake cache script *only if you intend to reproduce the CI build*.

### Additional Information

If you are building ExaGO on one of the following systems, you may follow the
link for more machine-specific information:

- [Ascent, Newell, or Marianas](./docs/web/README.ci_clusters.md)
- [Summit](./docs/web/README.summit.md)

