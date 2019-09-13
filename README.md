# SCOPFLOW
SCOPFLOW is an application code for security-constrained optimal power flow. It, currently, solves a two-stage stochastic optimization problem. The problem formulation is as follows:
```math
\text{min}&~\sum_{i=0}^{Ns} f(x_i)& \\
&\text{s.t.}& \\
&~g(x_i) = 0~~~i \in \{0,N_s\}& \\
&~h(x_i) \le 0~~i \in \{0,N_s\}& \\
-\delta{x}& \le x_i - x_0 \le \delta{x}~~i \in \{1,N_s\}&\\
&x^- \le x_i \le x^+~~i\in \{0,N_s\}&
 ```
where $N_s$ is the number of scenarios. Note that the last equation is the coupling for the 2nd stage scenarios and the first-stage.


## 1. Installing dependencies
SCOPFLOW is dependent on PETSc and either PIPS-NLP or IPOPT. Below, we provide installation instructions for these libraries.

### PETSc

#### Download
```
git clone https://bitbucket.org/petsc/petsc petsc
```
#### Set environment variables PETSC_ARCH and PETSC_DIR
PETSc requires two environment variables to be set to know the location (PETSC_DIR) and the configuration environment (PETSC_ARCH)
```
export PETSC_DIR=<petsc-location>
export PETSC_ARCH=<arch-name>
```
arch-name can be any name.

#### Installation
```
cd $PETSC_DIR
./config/configure.py --download-mpich --with-cc=gcc --with-\
cxx=g++ --with-fc=gfortran --download-mumps --download-scalapack --download-superlu --download-superlu_dist --download-suitespar\
se --download-metis --download-parmetis --download-cmake --with-cxx-dialect=C++11
make
make test
```
Run the above commands to install PETSc. Note that the installation options following configure.py are optional and can be omitted.

### PIPS-NLP
PIPS-NLP is a 2-stage parallel nonlinear stochastic optimization solver package. It uses an interior point method for the solution process and uses a Schur-complement approach to accelerate the solution. It is written in C++ with a C interface, which SCOPFLOW uses, available.

#### Download
```
git clone https://github.com/Argonne-National-Laboratory/PIPS.git
```

#### Install
See the installation instructions file README.md to install PIPS and its dependencies. Note that for SCOPFLOW, only PIPS-NLP will be used. Modify build_pips.sh file to disable building all libraries (-DBUILD_ALL=OFF) and enable PIPS-NLP (-DBUILD_PIPS_NLP)

#### Set environment variables
SCOPFLOW needs the location of the PIPS directory to include its header files and link with its libraries. This is done by setting an environment variables, PIPS_DIR, that points to the top-level PIPS directory location.
```
export PIPS_DIR=<path-to-pips>
```

### IPOPT
IPOPT is a popular nonlinear interior point method package for solving general nonlinear optimization problems. Like PIPS-NLP, it is written in C++ and has a C interface available.

#### Download
Visit the IPOPT site `https://github.com/coin-or/Ipopt` for download instructions for IPOPT. Note that IPOPT has several dependencies (in particular the HSL libraries) that need to be downloaded and put in the correct IPOPT folders before installing IPOPT. IPOPT tarballs are available at `https://www.coin-or.org/download/source/Ipopt/`


#### Install
Refer to the IPOPT installation instructions at `https://coin-or.github.io/Ipopt/INSTALL.html`

#### Set environment variables
SCOPFLOW needs to know the location of the IPOPT directory to correctly include its header files and link with the libraries. Set the environment variable IPOPT_BUILD_DIR to point to the location of IPOPT's build directory. Note that this the location of IPOPT's source code directory, it is the location where IPOPT was built.
```
export IPOPT_BUILD_DIR=<location-of-IPOPT-build-dir>
```

### Set environment variable for SCOPFLOW

SCOPFLOW needs the environment variable SCOPFLOW_DIR set to the location of the SCOPFLOW directory
```
export SCOPFLOW_DIR=<location of SCOPFLOW>
```

## 2. Compiling and running SCOPFLOW

### Compiling
SCOPFLOW currently uses separate codes for IPOPT and PIPS-NLP. In the future, we plan to unify these codes to share the common pieces. For now, there are separate executables for SCOPFLOW-IPOPT and SCOPFLOW-PIPSNLP implementation. To build SCOPFLOW with IPOPT do
```
make SCOPFLOW_IPOPT
```
and with PIPS-NLP do
```
make SCOPFLOW_PIPS
```
This will compile the main driver application for SCOPFLOW (applications/scopflow-main.c)
Note that if you have only want to use SCOPFLOW with IPOPT and have not installed PIPS-NLP then an additional compile option needs to set
```
make SCOPFLOW_IPOPT -DCFLAGS_PIPS=
```
and similiarly with PIPS-NLP it is
```
make SCOPFLOW_PIPS -DCFLAGS_IPOPT=
```

### Execution
Running SCOPFLOW-IPOPT
```
./SCOPFLOW_IPOPT
```
or SCOPFLOW-PIPS
```
mpiexec -n <N> ./SCOPFLOW_PIPS
```
## Options
The current version has several options available for SCOPFLOW. These options can be either set through the options file `scopflowoptions` or via the command line. Solver options are set through the files `ipopt.opt` and `pipsnlp.parameter` for IPOPT and PIPS-NLP, respectively.

#### Network file (-netfile <netfilename>): 
Set the name of the network file (MATPOWER format only currently). There is support for reading the PSSE format as well, but the PSSE raw data file does not contain the generator cost data. A separate file needs to be set for the generator cost which is not supported yet with SCOPFLOW.
```
./SCOPFLOW_IPOPT -netfile <netfilename>
```

#### Contingency file (-ctgcfile <ctgcfilename>): 
Set the name of the contingency data file. The contingency scenarios are set through this file. If this file is not set, then the base case scenario is repeated for each scenario.
```
./SCOPFLOW_IPOPT -netfile datafiles/case9.mod -ctgcfile <ctgcfilename>
```
The format of the contingencies is explained in the header file `include/scopflow.h`

#### Number of scenarios (-scopflow_Ns <Ns>): 
Sets the number of scenarios. This should be less than or equal to the number of contingencies set in the contingency file.
```
./SCOPFLOW_IPOPT -netfile datafiles/case9.mod -ctgcfile datafiles/case9.cont -scopflow_Ns <Ns>
```
With this option set, SCOPFLOW will only pick up the first Ns contingencies in the contingency file. 
With PIPS implementation, the number of scenarios set should be larger than the number of ranks.
```
mpiexec -n 3 ./SCOPFLOW_PIPS -netfile datafiles/case9.mod -ctgcfile datafiles/case9.cont -scopflow_Ns <Ns>
```
#### Coupling between first and second stage (-scopflow_iscoupling <0 or 1>)
```
./SCOPFLOW_IPOPT -netfile datafiles/case9.mod -ctgcfile datafiles/case9.cont -scopflow_coupling 1
```
#### Include first stage generation cost only (-scopflow_first_stage_gen_cost_only <0 or 1>)
```
./SCOPFLOW_IPOPT -netfile datafiles/case9.mod -ctgcfile datafiles/case9.cont -scopflow_first_stage_gen_cost_only 0
```

#### Ignore line flow constraints (-scopflow_ignore_line_flow_constraints <0 or 1>)
```
./SCOPFLOW_IPOPT -netfile datafiles/case9.mod -ctgcfile datafiles/case9.cont -scopflow_ignore_line_flow_constraints 1
```


### Additional installation tips
#### Linking
Depending on how you installed PIPS and IPOPT, you may need to update the linker environment variables to include the PETSc, IPOPT, PIPS-NLP, and SCOPFLOW library paths
```
export LD_LIBRARY_PATH=$PETSC_DIR/$PETSC_ARCH/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$SCOPFLOW_DIR:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$IPOPT_BUILD_DIR/lib:$LD_LIBRARY_PATH
LD_LIBRARY_PATH=$PIPS_DIR/build/PIPS-NLP:$LD_LIBRARY_PATH
```