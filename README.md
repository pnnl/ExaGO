# SCOPFLOW (<b>S</b>ecurity <b>C</b>onstrained <b>O</b>ptimal <b>P</b>ower <b>FLOW</b>)
SCOPFLOW is an application code for solving power grid security-constrained optimal power flow problem. The problem is set up as a two-stage optimization problem where the first-stage (base-case) represents the normal operation of the grid and the second-stage comprises of $N_s$ scenarios representing deviations from the normal operation such as those caused by contingencies or renewable uncertainties. Compactly, the problem can be set up in the following form:

```math
\begin{aligned}
\text{min}&~\sum_{i=0}^{Ns} f(x_i)& \\
&\text{s.t.}& \\
&~g(x_i) = 0~~~i \in \{0,N_s\}& \\
&~h(x_i) \le 0~~i \in \{0,N_s\}& \\
&x^- \le x_i \le x^+~~i\in \{0,N_s\}& \\
-\delta{x}& \le x_i - x_0 \le \delta{x}~~i \in \{1,N_s\}&
\end{aligned}
 ```
where $N_s$ is the number of scenarios. The last equation is the coupling between the 2nd stage scenarios and the first-stage that enforces a limit on the deviation on second stage decision variable $x_i$ from its corresponding base case decision variable $x_0$.

## 1. Download
SCOPFLOW can be downloaded from gitlab via
```
git clone https://gitlab.pnnl.gov/exasgd/frameworks/scopflow.git
```
or if you have SSH access to the repository
```
git clone ssh://git@gitlab.pnnl.gov:2222/exasgd/frameworks/scopflow.git
```

## 2. Install dependencies
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

### PIPS-NLP
PIPS-NLP is a 2-stage parallel nonlinear stochastic optimization solver package. It uses an interior point method for the solution process and uses a Schur-complement approach to accelerate the solution. It is written in C++ with a C interface, which SCOPFLOW uses, available. Note that currently the SCOPFLOW implementation with PIPS-NLP also depends on IPOPT so IPOPT needs to be also built using the above instructions. In the future, we plan to remove this dependency for SCOPFLOW PIPS-NLP implementation.

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

### Set environment variable for SCOPFLOW

SCOPFLOW needs the environment variable SCOPFLOW_DIR set to the location of the SCOPFLOW directory
```
export SCOPFLOW_DIR=<location of SCOPFLOW>
```

## 3. Building, execution, and options

### Compiling
SCOPFLOW can be either built with IPOPT or PIPS-NLP. To build with IPOPT do
```
make SCOPFLOW WITH_IPOPT=1
```
and for PIPS-NLP
```
make SCOPFLOW WITH_PIPS=1
```
This will compile the main driver application for SCOPFLOW (applications/scopflow-main.c)

### Execution
The SCOPFLOW code is executed via
```
mpiexec -n <N> ./SCOPFLOW <options>
```
where \<options\> are the available command line options as given in the next section.

### Options
The current version has several options available for SCOPFLOW. These options can be either set through the options file `options/scopflowoptions` or via the command line.

#### Network file (-netfile <netfilename>): 
Set the name of the network file (MATPOWER format only currently). There is support for reading the PSSE format as well, but the PSSE raw data file does not contain the generator cost data. A separate file needs to be set for the generator cost which is not supported yet with SCOPFLOW.

```
mpiexec -n <N> ./SCOPFLOW -netfile <netfilename>
```

#### Contingency file (-ctgcfile <ctgcfilename>): 
Set the name of the contingency data file. The contingency scenarios are set through this file.
```
mpiexec -n <N> ./SCOPFLOW -netfile <netfilename> -ctgcfile <ctgcfilename>
```
The format of the contingencies is explained in the header file `include/scopflow.h`

#### Solver (-scopflow_solver <IPOPT or PIPS>)
Set the solver to be used for SCOPFLOW. Currently, only IPOPT and PIPS-NLP are supported. If IPOPT is chosen then SCOPFLOW can be only run on one processor (N = 1) as IPOPT only supports single process execution. PIPS supports parallel execution.
```
mpiexec -n <N> ./SCOPFLOW -netfile <netfilename> -ctgcfile <ctgcfilename> -scopflow_solver <IPOPT or PIPS>
```
If this option is not set then SCOPFLOW uses IPOPT solver.

#### Formulation (-scopflow_formulation <formulationname>)
Set the formulation (representation of variables and equations) to be used for SCOPFLOW. The default formulation is power balance form with polar representation for voltages
```
mpiexec -n <N> ./SCOPFLOW -netfile <netfilename> -ctgcfile <ctgcfilename> -scopflow_formulation <formulationname>
```
Currently, three formulations are supported by SCOPFLOW.
1. POWER_BALANCE_POLAR: Power balance form with polar representation of voltages.
1. POWER_BALANCE_CARTESIAN: Power balance form with cartesian representation of voltages.
1. CURRENT_BALANCE_CARTESIAN: Current balance form with cartesian representation of voltages.

If this option is not set then SCOPFLOW uses POWER_BALANCE_POLAR formulation

#### Number of scenarios (-scopflow_Ns <Ns>): 
Sets the number of second-stage scenarios. This should be less than or equal to the number of contingencies set in the contingency file.
```
./SCOPFLOW -netfile <netfilename> -ctgcfile <ctgcfilename> -scopflow_Ns <Ns>
```
With this option set, SCOPFLOW will only pick up the first Ns contingencies in the contingency file. 
With PIPS as the solver for SCOPFLOW, the number of scenarios set should be larger than the number of ranks ($N_s > N$).

### Additional installation tips
#### Linking
Depending on how you installed PIPS and IPOPT, you may need to update the linker environment variables to include the PETSc, IPOPT, PIPS-NLP, and SCOPFLOW library paths
```
export LD_LIBRARY_PATH=$PETSC_DIR/$PETSC_ARCH/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$SCOPFLOW_DIR:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$IPOPT_BUILD_DIR/lib:$LD_LIBRARY_PATH
LD_LIBRARY_PATH=$PIPS_DIR/build/PIPS-NLP:$LD_LIBRARY_PATH
```