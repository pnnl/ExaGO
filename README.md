# SCOPFLOW
SCOPFLOW is an application code for security-constrained optimal power flow.

## Installation
SCOPFLOW makes heavy use of the PETSc library and its functionality.

### Download PETSc
```
git clone -b maint https://bitbucket.org/petsc/petsc petsc
```
### Set environment variables PETSC_ARCH and PETSC_DIR
PETSc requires two environment variables to be set to know the location (PETSC_DIR) and the configuration environment (PETSC_ARCH)
```
export PETSC_DIR=<petsc-location>
export PETSC_ARCH=<arch-name>
```
arch-name can be any name.

### Install PETSc
```
cd $PETSC_DIR
./config/configure.py --download-mpich --with-cc=gcc --with-\
cxx=g++ --with-fc=gfortran --download-mumps --download-scalapack --download-superlu --download-superlu_dist --download-suitespar\
se --download-metis --download-parmetis --download-cmake --with-cxx-dialect=C++11
make
make test
```
Run the above commands to install PETSc.

### Set environment variable for SCOPFLOW
SCOPFLOW needs the environment variable SCOPFLOW_DIR set to the location of the SCOPFLOW directory
```
export SCOPFLOW_DIR=<location of SCOPFLOW>
export LD_LIBRARY_PATH=$PETSC_DIR/$PETSC_ARCH/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$SCOPFLOW_DIR:$LD_LIBRARY_PATH

```
