/*
 * Public header file for security constrained optimal power flow application.
 *
 */

#ifndef SCOPFLOW_H
#define SCOPFLOW_H

#include <ps.h>

/* Models */
#define SCOPFLOWMODEL_GENRAMP "GENRAMP" 
#define SCOPFLOWMODEL_GENRAMPT "GENRAMPT" /* Model for multi-period SCOPFLOW */

/* Solvers */
#define SCOPFLOWSOLVER_IPOPTNEW "IPOPTNEW"
#define SCOPFLOWSOLVER_IPOPT "IPOPT"
#define SCOPFLOWSOLVER_TAO   "TAO"
#define SCOPFLOWSOLVER_EMPAR  "EMPAR" /* Embarassingly parallel solver - solves each OPFLOW independently */

typedef struct _p_SCOPFLOW *SCOPFLOW;

PETSC_EXTERN PetscErrorCode SCOPFLOWEnableMultiPeriod(SCOPFLOW,PetscBool);
PETSC_EXTERN PetscErrorCode SCOPFLOWSetModel(SCOPFLOW,const char[]);
PETSC_EXTERN PetscErrorCode SCOPFLOWModelRegister(SCOPFLOW,const char[],PetscErrorCode (*create)(SCOPFLOW));
PETSC_EXTERN PetscErrorCode SCOPFLOWSetSolver(SCOPFLOW,const char[]);
PETSC_EXTERN PetscErrorCode SCOPFLOWCreate(MPI_Comm,SCOPFLOW*);
PETSC_EXTERN PetscErrorCode SCOPFLOWDestroy(SCOPFLOW*);
PETSC_EXTERN PetscErrorCode SCOPFLOWSetNetworkData(SCOPFLOW,const char[]);

/*
 * The native format for Contingency input file is as follows:
 * Notes: Each field in the contingency file has the following format
 * Num,Type,Bus,Fbus,Tbus,Id,Status,prob
 * Num - Contingency number
 * Type - Type of contingency (Generator, Branch, Transformer, Load)
 * Bus - The equipment bus number
 * Fbus - From bus number (only for branch and transformer contingencies)
 * Tbus - To bus number (only for branch and transformer contingencies)
 * Id   - The equipment ID (2-char string)
 * Status - The status to be set for the outaged equipment
 * Prob   - the probability of the outage

  Examples:
(Outage generator at bus 1 with Id 1, probability = 0.1)
    1,0,1,0,0,1 ,0,0.1  
(Outage branch connecting buses 8-9 with Id 1, probability = 0.1)
    2,1,0,8,9,1 ,0,0.1
(Multiple outages
  generator at bus 1 with Id 1, probability = 0.1
  branch connecting buses 8-9 with Id 1, probability = 0.1)
    3,0,1,0,0,1 ,0,0.1  
    3,1,0,8,9,1 ,0,0.1
 */
PETSC_EXTERN PetscErrorCode SCOPFLOWSetContingencyData(SCOPFLOW,ContingencyFileInputFormat,const char[]);
PETSC_EXTERN PetscErrorCode SCOPFLOWSetUp(SCOPFLOW);
PETSC_EXTERN PetscErrorCode SCOPFLOWCreateGlobalVector(SCOPFLOW,Vec*);
PETSC_EXTERN PetscErrorCode SCOPFLOWCreateMatrix(SCOPFLOW,Mat*);
PETSC_EXTERN PetscErrorCode SCOPFLOWSolve(SCOPFLOW);
PETSC_EXTERN PetscErrorCode SCOPFLOWSetNumScenarios(SCOPFLOW,PetscInt);
PETSC_EXTERN PetscErrorCode SCOPFLOWGetObjective(SCOPFLOW,PetscReal*);
PETSC_EXTERN PetscErrorCode SCOPFLOWGetSolution(SCOPFLOW,PetscInt,Vec*);
PETSC_EXTERN PetscErrorCode SCOPFLOWPrintSolution(SCOPFLOW,PetscInt);
PETSC_EXTERN PetscErrorCode SCOPFLOWSaveSolution(SCOPFLOW,PetscInt,OutputFormat,const char[]);
PETSC_EXTERN PetscErrorCode SCOPFLOWSaveSolutionAll(SCOPFLOW,OutputFormat,const char[]);
PETSC_EXTERN PetscErrorCode SCOPFLOWGetConvergenceStatus(SCOPFLOW,PetscBool*);
PETSC_EXTERN PetscErrorCode SCOPFLOWGetMode(SCOPFLOW,PetscInt*);

#endif


