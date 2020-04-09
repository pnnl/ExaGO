## PIPS-NLP

PIPS-NLP is a 2-stage parallel nonlinear stochastic optimization solver package. It uses an interior point method for the solution process and uses a Schur-complement approach to accelerate the solution. It is written in C++ with a C interface, which SCOPFLOW uses, available. Note that currently the SCOPFLOW implementation with PIPS-NLP also depends on IPOPT so IPOPT needs to be also built using the above instructions. In the future, we plan to remove this dependency for SCOPFLOW PIPS-NLP implementation.

### Download
```
git clone https://github.com/Argonne-National-Laboratory/PIPS.git
```

### Install
See the installation instructions file README.md to install PIPS and its dependencies. Note that for SCOPFLOW, only PIPS-NLP will be used. Modify build_pips.sh file to disable building all libraries (-DBUILD_ALL=OFF) and enable PIPS-NLP (-DBUILD_PIPS_NLP=ON).
PIPS-NLP also has a dependency on the partitioning libraries METIS and ParMETIS. Since these have already been installed with PETSc, they can be reused with PIPS-NLP. Modify build_pips.sh to add a flag
(-DMETIS_DIR=$PETSC_DIR/$PETSC_ARCH) to use METIS/ParMETIS installed with PETSc. We recommend using METIS/ParMetis installed with PETSc instead of
using the METIS downloaded via PIPS as this will avoid some of the library conflicts that may happen. SCOPFLOW assumes that PIPS is installed in $PIPS_DIR/build_pips (as given in the PIPS README.md file)
### Set environment variables (only with `make` build)
SCOPFLOW needs the location of the PIPS directory to include its header files and link with its libraries. This is done by setting an environment variables, PIPS_DIR, that points to the top-level PIPS directory location.
```
export PIPS_DIR=<path-to-pips>
```