## PETSc

These instructions need to be updated to expand on the myriad of ways PETSc can be built.

### Download
Download the PETSc release version (ver 3.13) 
```
git clone https://bitbucket.org/petsc/petsc petsc
```
### Set environment variables PETSC_ARCH and PETSC_DIR
PETSc requires two environment variables to be set to know the location (PETSC_DIR) and the configuration environment (PETSC_ARCH)
```
export PETSC_DIR=<petsc-location>
export PETSC_ARCH=<arch-name>
```
arch-name can be any name.

### Switch branch
SCOPFLOW is compatible with the current release of PETSc (version 3.13). To use it, switch to the `maint` branch

### Installation
```
cd $PETSC_DIR
git checkout maint
./config/configure.py --download-mpich --with-cc=gcc --with-\
cxx=g++ --with-fc=gfortran --download-mumps --download-scalapack --download-superlu --download-superlu_dist --download-suitespar\
se --download-metis --download-parmetis --download-cmake --with-cxx-dialect=C++11
make
make check
```
Run the above commands to install PETSc. Note that the installation options following configure.py are optional and can be omitted.