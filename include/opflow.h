/**
 * @file opflow.h
 * @brief Public header file for optimal power flow application.
 *
 */

#ifndef OPFLOW_H
#define OPFLOW_H

#include <ps.h>

typedef struct _p_OPFLOW *OPFLOW;

PETSC_EXTERN PetscErrorCode OPFLOWSetFormulation(OPFLOW,const char*);
PETSC_EXTERN PetscErrorCode OPFLOWSetSolver(OPFLOW,const char*);
PETSC_EXTERN PetscErrorCode OPFLOWFormulationRegister(const char[],PetscErrorCode (*create)(OPFLOW));
PETSC_EXTERN PetscErrorCode OPFLOWSolverRegister(const char[],PetscErrorCode (*create)(OPFLOW));
PETSC_EXTERN PetscErrorCode OPFLOWCreate(MPI_Comm,OPFLOW*);
PETSC_EXTERN PetscErrorCode OPFLOWDestroy(OPFLOW*);
PETSC_EXTERN PetscErrorCode OPFLOWReadMatPowerData(OPFLOW,const char[]);
PETSC_EXTERN PetscErrorCode OPFLOWSetUp(OPFLOW);
PETSC_EXTERN PetscErrorCode OPFLOWCreateGlobalVector(OPFLOW,Vec*);
PETSC_EXTERN PetscErrorCode OPFLOWCreateMatrix(OPFLOW,Mat*);
PETSC_EXTERN PetscErrorCode OPFLOWSolve(OPFLOW);
#endif


