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
#define OPFLOWMODEL_PBPOLHIOP "POWER_BALANCE_HIOP"
#define OPFLOWMODEL_PBPOLRAJAHIOP "PBPOLRAJAHIOP"

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
PETSC_EXTERN PetscErrorCode OPFLOWGetConstraintJacobian(OPFLOW,Mat*,Mat*);
PETSC_EXTERN PetscErrorCode OPFLOWGetHessian(OPFLOW,Mat*,PetscScalar*);
PETSC_EXTERN PetscErrorCode OPFLOWGetSolution(OPFLOW,Vec*);
PETSC_EXTERN PetscErrorCode OPFLOWGetConvergenceStatus(OPFLOW,PetscBool*);
PETSC_EXTERN PetscErrorCode OPFLOWGetConstraints(OPFLOW,Vec*);
PETSC_EXTERN PetscErrorCode OPFLOWGetConstraintMultipliers(OPFLOW,Vec*);

PETSC_EXTERN PetscErrorCode OPFLOWComputeVariableBounds(OPFLOW,Vec,Vec);
PETSC_EXTERN PetscErrorCode OPFLOWComputeGradient(OPFLOW,Vec,Vec);
PETSC_EXTERN PetscErrorCode OPFLOWComputeGradientArray(OPFLOW,const double*,double*);
PETSC_EXTERN PetscErrorCode OPFLOWComputeObjective(OPFLOW,Vec,PetscReal*);
PETSC_EXTERN PetscErrorCode OPFLOWComputeObjectiveArray(OPFLOW,const double*,double*);
PETSC_EXTERN PetscErrorCode OPFLOWComputeConstraints(OPFLOW,Vec,Vec);
PETSC_EXTERN PetscErrorCode OPFLOWComputeEqualityConstraints(OPFLOW,Vec,Vec);
PETSC_EXTERN PetscErrorCode OPFLOWComputeInequalityConstraints(OPFLOW,Vec,Vec);
PETSC_EXTERN PetscErrorCode OPFLOWComputeEqualityConstraintsArray(OPFLOW,const double*,double*);
PETSC_EXTERN PetscErrorCode OPFLOWComputeInequalityConstraintsArray(OPFLOW,const double*,double*);

PETSC_EXTERN PetscErrorCode OPFLOWGetConstraintBounds(OPFLOW,Vec*,Vec*);
PETSC_EXTERN PetscErrorCode OPFLOWComputeConstraintBounds(OPFLOW,Vec,Vec);
PETSC_EXTERN PetscErrorCode OPFLOWComputeConstraintJacobian(OPFLOW,Vec,Mat,Mat);
PETSC_EXTERN PetscErrorCode OPFLOWComputeHessian(OPFLOW,Vec,Vec,PetscScalar,Mat);

PETSC_EXTERN PetscErrorCode OPFLOWPrintSolution(OPFLOW);
PETSC_EXTERN PetscErrorCode OPFLOWSaveSolution(OPFLOW,OutputFormat,const char[]);
PETSC_EXTERN PetscErrorCode OPFLOWGenbusVoltageFixed(OPFLOW,PetscBool);

PETSC_EXTERN PetscErrorCode OPFLOWGetVariableOrdering(OPFLOW,int**);
PETSC_EXTERN PetscErrorCode OPFLOWGetSizes(OPFLOW,int*,int*,int*);
#endif


