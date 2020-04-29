## PETSc

Portable Extensible Toolkit for Scientific Computation ([PETSc](https://www.mcs.anl.gov/petsc/)) is a popular high-performance numerical library for large-scale scientific computation. It includes a suite of high-performance linear, nonlinear, time-stepping, and optimization solvers. PETSc also includes a framework for managing unstructured networks called [`DMNetwork`](https://www.mcs.anl.gov/petsc/dmnetwork/index.html) which is used by SCOPFLOW for topology and data (component) management.

### 1. Download
Download the PETSc release version (ver 3.13) 
```
git clone https://bitbucket.org/petsc/petsc petsc
```
### 2. Set environment variables PETSC_ARCH and PETSC_DIR
PETSc requires two environment variables to be set to know the location (PETSC_DIR) and the configuration environment (PETSC_ARCH)
```
export PETSC_DIR=<petsc-location>
export PETSC_ARCH=<arch-name>
```
arch-name can be any name.

### 3. Switch branch
SCOPFLOW is compatible with the current release of PETSc (version 3.13). To use it, switch to the `maint` branch

```
cd $PETSC_DIR
git checkout maint
```

### 4. Configuration
PETSc can be installed in myriad of ways. In-depth details on PETSc installation are given [here](https://www.mcs.anl.gov/petsc/documentation/installation.html). Below, we describe a few ways to install PETSc and (optionally) install different packages with it.
#### Minimal
In its minimal configuration, only the configure script needs to be run.
```
./config/configure.py
```
In this mode, PETSc is not installed with any third party packages. It will look for MPI and BLAS/LAPACK (its core dependencies) in certain standard locations. Since SCOPFLOW has only optional dependence on third party packages installed with PETSc, running just the configure script should suffice.

By default, PETSc is installed in debug mode. To use the release/production mode, add the flag `--with-debugging=0`. Compiler optimization flags can also be specified if need be.

```
./config/configure.py --with-debugging=0 COPTFLAGS=-O3 CXXOPTFLAGS=-O3
```

To change the location of the PETSc installation use `--prefix=<location_of_PETSC_install>`. If prefix option is not set then PETSc will install files at `/usr`. We highly recommend setting a PETSc installation directory with the `--prefix` option.

#### With MPI and/or BLAS/LAPACK specified
PETSc will attempt to find its core dependencies (MPI, Blas/Lapack) in certain standard locations. If these are in non-standard locations then they can be specified via 
- MPI: `--with-mpi-dir=<location_of_mpi>`
- BLAS/LAPACK: `--with-blaslapack-dir=<location_of_blaslapack>`

Alternately, one can also download MPI and/or BLAS/LAPACK
- MPI: `--download-mpich`
- BLAS/LAPACK: `--download-fblaslapack` (when fortran compiler is present) or `--download-f2cblaslapack` (without a fortran compiler)

#### With third-party packages
PETSc has interfaces to a variety of third party packages. These third party packages can be either downloaded and installed with PETSc installation, or if already installed can be linked to PETSc. We recommend downloading the packages to avoid any conflict or version incompatability issues. Below, we specify a few of these solver packages

- METIS:
    - Download: `--download-metis`
    - Already installed: `--with-metis-dir=<metis_location>`
- PARMETIS:
    - Download: `--download-parmetis`
    - Already installed: `--with-parmetis-dir=<metis_location>`

If you are installing PARMETIS then installing METIS is also suggested in case there are any dependency issues.

- SuperLU\_Dist:
    - Download: `--download-superlu_dist --download-metis --download-parmetis`
    - Already installed: `--with-superlu_dist-dir=<superlu_dist_location> --with-metis-dir=<metis_location> --with-parmetis-dir=<parmetis_location>` 
    
    Note that SuperLU_Dist depends on METIS and PARMETIS libraries hence these are needed as well.

- MUMPS:
    - Download: `--download-mumps --download-metis --download-parmetis --download-scalapack`
    - Already installed: `--with-mumps-dir=<superlu_dist_location> --with-metis-dir=<metis_location> --with-parmetis-dir=<parmetis_location> --with-scalapack-dir=<scalapack_location>` 
    
    Note that SuperLU_Dist depends on METIS and PARMETIS libraries hence these are needed as well.

- SuiteSparse:
    - Download: `--download-suitesparse`
    - Already installed: `--with-suitesparse-dir=<suitesparse_location>`

    Note that PETSc can only use AMD, CHOLMOD, UMFPACK, and KLU from the SuiteSparse collection.

#### Example full install

```
./config/configure.py --download-mpich --download-mumps --download-scalapack --download-superlu_dist --download-suitesparse --download-metis --download-parmetis
```

### 5. Compile, check, and install
Once configure goes through successfully, compile and check PETSc installation with
```
make
make check
make install
```
At the completion of this step, PETSc installation files will be at the location set with `--prefix` configure option. If `--prefix` option is not set then PETSc will be installed in `/usr`.

