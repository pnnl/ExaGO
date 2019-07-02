# PSapps
PSapps is a library of high-performance power grid applications. It provides two applications currently
- Power Flow (PFLOW)
- Dynamics Simulation (DYN)

## Installation
PSAPPS makes heavy use of the PETSc library and its functionality.

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

### Set environment variable for PSapps
PSapps needs the environment variable PSAPPS_DIR set to the location of the PSapps directory
```
export PSAPPS_DIR=<location of PSapps>
```
## Power flow application
Compile the power flow application code
```
make PFLOW
```
This will create the executable PFLOW in the PSAPPS_DIR. The source code for the PFLOW application is in the `applications` directory pflow-main.c

```
mpiexec -n <n> ./PFLOW -netfile <networkfile_name>
```
Here <networkfile_name> is the name of network file. PFLOW can read power grid networks in MATPOWER and PSSE RAW data formats.
All PETSc solver options can be passed in with the prefix -pflow_. For example the PETSc option -pc_type would -pflow_pc_type

PFLOW reads the options file pflowoptions located in PSAPPS_DIR. All options can be set in this file.
