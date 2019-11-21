/**
 * @file scopflow.h
 * @brief Public header file for security constrained optimal power flow application.
 *
 */

#ifndef SCOPFLOW_H
#define SCOPFLOW_H

#include <ps.h>

/* Formulations */
#define SCOPFLOWFORMULATION_PBPOL "POWER_BALANCE_POLAR" 
#define SCOPFLOWFORMULATION_PBCAR "POWER_BALANCE_CARTESIAN"
#define SCOPFLOWFORMULATION_IBCAR "CURRENT_BALANCE_CARTESIAN"

/* Solvers */
#define SCOPFLOWSOLVER_IPOPT "IPOPT"
#define SCOPFLOWSOLVER_TAO   "TAO"
#define SCOPFLOWSOLVER_PIPS  "PIPS"


typedef struct _p_SCOPFLOW *SCOPFLOW;

PETSC_EXTERN PetscErrorCode SCOPFLOWSetFormulation(SCOPFLOW,const char[]);
PETSC_EXTERN PetscErrorCode SCOPFLOWSetSolve(SCOPFLOW,const char[]);
PETSC_EXTERN PetscErrorCode SCOPFLOWCreate(MPI_Comm,SCOPFLOW*);
PETSC_EXTERN PetscErrorCode SCOPFLOWDestroy(SCOPFLOW*);
PETSC_EXTERN PetscErrorCode SCOPFLOWSetNetworkData(SCOPFLOW,const char[]);

/**
 * @brief Sets the contingency data file 
 * @param [in] SCOPFLOW scopflow - The SCOPFLOW object
 * @param [in] const char[] ctgcfile - The name of the contingency list file
 *
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
PETSC_EXTERN PetscErrorCode SCOPFLOWSetContingencyData(SCOPFLOW,const char[]);
PETSC_EXTERN PetscErrorCode SCOPFLOWSetUp(SCOPFLOW);
PETSC_EXTERN PetscErrorCode SCOPFLOWCreateGlobalVector(SCOPFLOW,Vec*);
PETSC_EXTERN PetscErrorCode SCOPFLOWCreateMatrix(SCOPFLOW,Mat*);
PETSC_EXTERN PetscErrorCode SCOPFLOWSolve(SCOPFLOW);
PETSC_EXTERN PetscErrorCode SCOPFLOWSetNumScenarios(SCOPFLOW,PetscInt);

#endif


