/*
 * Public header file for optimal power flow application.
 */

#ifndef OPFLOW_H
#define OPFLOW_H

#include <ps.h>

/* Models */
#define OPFLOWMODEL_PBPOL "POWER_BALANCE_POLAR" 
#define OPFLOWMODEL_PBCAR "POWER_BALANCE_CARTESIAN"
#define OPFLOWMODEL_IBCAR "CURRENT_BALANCE_CARTESIAN"
#define OPFLOWMODEL_IBCAR2 "CURRENT_BALANCE_CARTESIAN2"
#define OPFLOWMODEL_PBPOL2 "POWER_BALANCE_POLAR2"

/* Solvers */
#define OPFLOWSOLVER_IPOPT "IPOPT"
#define OPFLOWSOLVER_TAO   "TAO"
#define OPFLOWSOLVER_HIOP  "HIOP"

typedef struct _p_OPFLOW *OPFLOW;

PETSC_EXTERN PetscErrorCode OPFLOWSetModel(OPFLOW,const char[]);
PETSC_EXTERN PetscErrorCode OPFLOWSetSolver(OPFLOW,const char[]);
PETSC_EXTERN PetscErrorCode OPFLOWModelRegister(OPFLOW,const char[],PetscErrorCode (*create)(OPFLOW));
PETSC_EXTERN PetscErrorCode OPFLOWSolverRegister(OPFLOW,const char[],PetscErrorCode (*create)(OPFLOW));
PETSC_EXTERN PetscErrorCode OPFLOWCreate(MPI_Comm,OPFLOW*);
PETSC_EXTERN PetscErrorCode OPFLOWDestroy(OPFLOW*);
PETSC_EXTERN PetscErrorCode OPFLOWReadMatPowerData(OPFLOW,const char[]);
PETSC_EXTERN PetscErrorCode OPFLOWSetUp(OPFLOW);
PETSC_EXTERN PetscErrorCode OPFLOWCreateGlobalVector(OPFLOW,Vec*);
PETSC_EXTERN PetscErrorCode OPFLOWCreateMatrix(OPFLOW,Mat*);
PETSC_EXTERN PetscErrorCode OPFLOWSolve(OPFLOW);
PETSC_EXTERN PetscErrorCode OPFLOWSetInitialGuess(OPFLOW,Vec);
PETSC_EXTERN PetscErrorCode OPFLOWGetObjective(OPFLOW,PetscReal*);
PETSC_EXTERN PetscErrorCode OPFLOWGetVariableBounds(OPFLOW,Vec*,Vec*);
PETSC_EXTERN PetscErrorCode OPFLOWComputeVariableBounds(OPFLOW,Vec,Vec);
PETSC_EXTERN PetscErrorCode OPFLOWComputeGradient(OPFLOW,Vec,Vec);
PETSC_EXTERN PetscErrorCode OPFLOWComputeObjective(OPFLOW,Vec,PetscReal*);
PETSC_EXTERN PetscErrorCode OPFLOWGetSolution(OPFLOW,Vec*);
PETSC_EXTERN PetscErrorCode OPFLOWGetConvergenceStatus(OPFLOW,PetscBool*);
PETSC_EXTERN PetscErrorCode OPFLOWGetConstraints(OPFLOW,Vec*);
PETSC_EXTERN PetscErrorCode OPFLOWComputeConstraints(OPFLOW,Vec,Vec);
PETSC_EXTERN PetscErrorCode OPFLOWGetConstraintBounds(OPFLOW,Vec*,Vec*);
PETSC_EXTERN PetscErrorCode OPFLOWComputeConstraintBounds(OPFLOW,Vec,Vec);
PETSC_EXTERN PetscErrorCode OPFLOWGetConstraintMultipliers(OPFLOW,Vec*);
PETSC_EXTERN PetscErrorCode OPFLOWPrintSolution(OPFLOW);
PETSC_EXTERN PetscErrorCode OPFLOWSaveSolution(OPFLOW,OutputFormat,const char[]);
#endif


