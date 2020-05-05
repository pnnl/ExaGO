/*
 * Public header file for security constrained optimal power flow application.
 *
 */

#ifndef TCOPFLOW_H
#define TCOPFLOW_H

#include <ps.h>

/* Solvers */
#define TCOPFLOWSOLVER_IPOPT "IPOPT"

typedef struct _p_TCOPFLOW *TCOPFLOW;

PETSC_EXTERN PetscErrorCode TCOPFLOWSetFormulation(TCOPFLOW,const char[]);
PETSC_EXTERN PetscErrorCode TCOPFLOWSetSolve(TCOPFLOW,const char[]);
PETSC_EXTERN PetscErrorCode TCOPFLOWCreate(MPI_Comm,TCOPFLOW*);
PETSC_EXTERN PetscErrorCode TCOPFLOWDestroy(TCOPFLOW*);
PETSC_EXTERN PetscErrorCode TCOPFLOWSetNetworkData(TCOPFLOW,const char[]);

PETSC_EXTERN PetscErrorCode TCOPFLOWSetUp(TCOPFLOW);
PETSC_EXTERN PetscErrorCode TCOPFLOWSolve(TCOPFLOW);
PETSC_EXTERN PetscErrorCode TCOPFLOWSetTimeStepandDuration(TCOPFLOW,PetscReal,PetscReal);
PETSC_EXTERN PetscErrorCode TCOPFLOWSetLoadProfiles(TCOPFLOW, const char[], const char[]);
PETSC_EXTERN PetscErrorCode TCOPFLOWSetWindGenProfiles(TCOPFLOW, const char[]);

#endif


