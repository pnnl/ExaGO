/*
 * Public header file for security constrained optimal power flow application.
 *
 */

#ifndef TCOPFLOW_H
#define TCOPFLOW_H

#include <ps.h>

/* Solvers */
#define TCOPFLOWSOLVER_IPOPT "IPOPT"

/* Models */
#define TCOPFLOWMODEL_GENRAMP "GENRAMP"

/* Initialization and Parameters*/
#define TCOPFLOW_INITIALIZATION "ACPF"
#define TCOPFLOW_GENBUSVOLTAGE "VARIABLE_WITHIN_BOUNDS"

typedef struct _p_TCOPFLOW *TCOPFLOW;

PETSC_EXTERN PetscErrorCode TCOPFLOWSetModel(TCOPFLOW,const char[]);
PETSC_EXTERN PetscErrorCode TCOPFLOWModelRegister(TCOPFLOW,const char[],PetscErrorCode (*create)(TCOPFLOW));
PETSC_EXTERN PetscErrorCode TCOPFLOWSetSolver(TCOPFLOW,const char[]);
PETSC_EXTERN PetscErrorCode TCOPFLOWCreate(MPI_Comm,TCOPFLOW*);
PETSC_EXTERN PetscErrorCode TCOPFLOWDestroy(TCOPFLOW*);
PETSC_EXTERN PetscErrorCode TCOPFLOWSetNetworkData(TCOPFLOW,const char[]);

PETSC_EXTERN PetscErrorCode TCOPFLOWSetUp(TCOPFLOW);
PETSC_EXTERN PetscErrorCode TCOPFLOWSolve(TCOPFLOW);
PETSC_EXTERN PetscErrorCode TCOPFLOWSetTimeStepandDuration(TCOPFLOW,PetscReal,PetscReal);
PETSC_EXTERN PetscErrorCode TCOPFLOWSetLoadProfiles(TCOPFLOW, const char[], const char[]);
PETSC_EXTERN PetscErrorCode TCOPFLOWSetWindGenProfiles(TCOPFLOW, const char[]);

PETSC_EXTERN PetscErrorCode TCOPFLOWGetObjective(TCOPFLOW,PetscReal*);
PETSC_EXTERN PetscErrorCode TCOPFLOWGetSolution(TCOPFLOW,PetscInt,Vec*);
PETSC_EXTERN PetscErrorCode TCOPFLOWPrintSolution(TCOPFLOW,PetscInt);
PETSC_EXTERN PetscErrorCode TCOPFLOWSaveSolution(TCOPFLOW,PetscInt,OutputFormat,const char[]);
PETSC_EXTERN PetscErrorCode TCOPFLOWSaveSolutionAll(TCOPFLOW,OutputFormat,const char[]);
PETSC_EXTERN PetscErrorCode TCOPFLOWGetConvergenceStatus(TCOPFLOW,PetscBool*);
PETSC_EXTERN PetscErrorCode TCOPFLOWGetNumIterations(TCOPFLOW,PetscInt*);
PETSC_EXTERN PetscErrorCode TCOPFLOWSetTolerance(TCOPFLOW,PetscReal);
PETSC_EXTERN PetscErrorCode TCOPFLOWGetTolerance(TCOPFLOW,PetscReal*);

#endif


