/*
 * Public header file for power flow application
 */

#ifndef PFLOW_H
#define PFLOW_H

#include <ps.h>
#include <utils.h>

typedef struct _p_PFLOW *PFLOW;

PETSC_EXTERN PetscErrorCode PFLOWCreate(MPI_Comm, PFLOW *);
PETSC_EXTERN PetscErrorCode PFLOWDestroy(PFLOW *);
PETSC_EXTERN PetscErrorCode PFLOWReadMatPowerData(PFLOW, const char[]);
PETSC_EXTERN PetscErrorCode PFLOWReadPSSERawData(PFLOW, const char[]);
PETSC_EXTERN PetscErrorCode PFLOWSetUp(PFLOW);
PETSC_EXTERN PetscErrorCode PFLOWCreateGlobalVector(PFLOW, Vec *);
PETSC_EXTERN PetscErrorCode PFLOWCreateMatrix(PFLOW, Mat *);
PETSC_EXTERN PetscErrorCode PFLOWSolve(PFLOW);
PETSC_EXTERN PetscErrorCode PFLOWSolutionToPS(PFLOW);
PETSC_EXTERN PetscErrorCode PFLOWConverged(PFLOW, PetscBool *);
PETSC_EXTERN PetscErrorCode PFLOWSetLineStatus(PFLOW, PetscInt, PetscInt,
                                               const char *, PetscInt);
PETSC_EXTERN PetscErrorCode PFLOWSetGenStatus(PFLOW, PetscInt, const char *,
                                              PetscInt);
PETSC_EXTERN PetscErrorCode PFLOWGetBusVoltage(PFLOW, PetscInt, PetscScalar *,
                                               PetscScalar *, PetscBool *);
PETSC_EXTERN PetscErrorCode PFLOWSetBusVoltage(PFLOW, PetscInt, PetscScalar,
                                               PetscScalar);
PETSC_EXTERN PetscErrorCode PFLOWGetGenDispatch(PFLOW, PetscInt, PetscInt,
                                                PetscScalar *, PetscScalar *,
                                                PetscBool *);
PETSC_EXTERN PetscErrorCode PFLOWSetGenDispatch(PFLOW, PetscInt, PetscInt,
                                                PetscScalar, PetscScalar);
PETSC_EXTERN PetscErrorCode PFLOWSetLoadPower(PFLOW, PetscInt, PetscScalar,
                                              PetscScalar);
PETSC_EXTERN PetscErrorCode PFLOWGetLoadPower(PFLOW, PetscInt, PetscScalar *,
                                              PetscScalar *, PetscBool *);
PETSC_EXTERN PetscErrorCode PFLOWAddBusShunt(PFLOW, PetscInt, PetscScalar,
                                             PetscScalar);
PETSC_EXTERN PetscErrorCode PFLOWSetInitialGuess(PFLOW, Vec);
PETSC_EXTERN PetscErrorCode PFLOWGetNumIterations(PFLOW, PetscInt *);
PETSC_EXTERN PetscErrorCode PFLOWGetConvergenceStatus(PFLOW,PetscBool *);
PETSC_EXTERN PetscErrorCode PFLOWPrintSolution(PFLOW);
PETSC_EXTERN PetscErrorCode PFLOWSetGICData(PFLOW, const char[]);
PETSC_EXTERN PetscErrorCode PFLOWSaveSolutionDefault(PFLOW, const char[]);
PETSC_EXTERN PetscErrorCode PFLOWSetOutputFormat(PFLOW, OutputFormat);

#endif
