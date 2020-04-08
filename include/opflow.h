/**
 * @file opflow.h
 * @brief Public header file for optimal power flow application.
 *
 */

#ifndef OPFLOW_H
#define OPFLOW_H

#include <ps.h>

/* Formulations */
#define OPFLOWFORMULATION_PBPOL "POWER_BALANCE_POLAR" 
#define OPFLOWFORMULATION_PBCAR "POWER_BALANCE_CARTESIAN"
#define OPFLOWFORMULATION_IBCAR "CURRENT_BALANCE_CARTESIAN"
#define OPFLOWFORMULATION_IBCAR2 "CURRENT_BALANCE_CARTESIAN2"

/* Solvers */
#define OPFLOWSOLVER_IPOPT "IPOPT"
#define OPFLOWSOLVER_TAO   "TAO"
#define OPFLOWSOLVER_HIOP  "HIOP"

typedef struct _p_OPFLOW *OPFLOW;

PETSC_EXTERN PetscErrorCode OPFLOWSetFormulation(OPFLOW,const char[]);
PETSC_EXTERN PetscErrorCode OPFLOWSetSolver(OPFLOW,const char[]);
PETSC_EXTERN PetscErrorCode OPFLOWFormulationRegister(OPFLOW,const char[],PetscErrorCode (*create)(OPFLOW));
PETSC_EXTERN PetscErrorCode OPFLOWSolverRegister(OPFLOW,const char[],PetscErrorCode (*create)(OPFLOW));
PETSC_EXTERN PetscErrorCode OPFLOWCreate(MPI_Comm,OPFLOW*);
PETSC_EXTERN PetscErrorCode OPFLOWDestroy(OPFLOW*);
PETSC_EXTERN PetscErrorCode OPFLOWReadMatPowerData(OPFLOW,const char[]);
PETSC_EXTERN PetscErrorCode OPFLOWSetUp(OPFLOW);
PETSC_EXTERN PetscErrorCode OPFLOWCreateGlobalVector(OPFLOW,Vec*);
PETSC_EXTERN PetscErrorCode OPFLOWCreateMatrix(OPFLOW,Mat*);
PETSC_EXTERN PetscErrorCode OPFLOWSolve(OPFLOW);
PETSC_EXTERN PetscErrorCode OPFLOWSetInitialGuess(OPFLOW,Vec);
#endif


