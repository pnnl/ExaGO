/*
 * Public header file for stochastic optimal power flow application.
 *
 */

#ifndef SOPFLOW_H
#define SOPFLOW_H

#include <ps.h>

/* Models */
#define SOPFLOWMODEL_GENRAMP "GENRAMP" 
#define SOPFLOWMODEL_GENRAMPC "GENRAMPC" /* Model for multi-contingency SOPFLOW */


/* Solvers */
#define SOPFLOWSOLVER_IPOPT "IPOPT"
#define SOPFLOWSOLVER_EMPAR  "EMPAR" /* Embarassingly parallel solver - solves each OPFLOW independently */
#define SOPFLOWSOLVER_HIOP   "HIOP" /* Primal decomposition-basedd HIOp solver */

/* Initialization and Parameters*/
#define SOPFLOW_INITIALIZATION "ACPF"
#define SOPFLOW_GENBUSVOLTAGE "VARIABLE_WITHIN_BOUNDS"

/* Type of uncertainty */
typedef enum {NONE,WIND,LOAD}ScenarioUncertaintyType;

typedef struct _p_SOPFLOW *SOPFLOW;

PETSC_EXTERN PetscErrorCode SOPFLOWEnableMultiContingency(SOPFLOW,PetscBool);
PETSC_EXTERN PetscErrorCode SOPFLOWSetModel(SOPFLOW,const char[]);
PETSC_EXTERN PetscErrorCode SOPFLOWModelRegister(SOPFLOW,const char[],PetscErrorCode (*create)(SOPFLOW));
PETSC_EXTERN PetscErrorCode SOPFLOWSetSolver(SOPFLOW,const char[]);
PETSC_EXTERN PetscErrorCode SOPFLOWCreate(MPI_Comm,SOPFLOW*);
PETSC_EXTERN PetscErrorCode SOPFLOWDestroy(SOPFLOW*);
PETSC_EXTERN PetscErrorCode SOPFLOWSetNetworkData(SOPFLOW,const char[]);

PETSC_EXTERN PetscErrorCode SOPFLOWSetScenarioData(SOPFLOW,ScenarioFileInputFormat,ScenarioUncertaintyType,const char[]);
PETSC_EXTERN PetscErrorCode SOPFLOWSetContingencyData(SOPFLOW,ContingencyFileInputFormat,const char[]);
PETSC_EXTERN PetscErrorCode SOPFLOWSetUp(SOPFLOW);
PETSC_EXTERN PetscErrorCode SOPFLOWCreateGlobalVector(SOPFLOW,Vec*);
PETSC_EXTERN PetscErrorCode SOPFLOWCreateMatrix(SOPFLOW,Mat*);
PETSC_EXTERN PetscErrorCode SOPFLOWSolve(SOPFLOW);
PETSC_EXTERN PetscErrorCode SOPFLOWSetNumScenarios(SOPFLOW,PetscInt);
PETSC_EXTERN PetscErrorCode SOPFLOWGetNumScenarios(SOPFLOW,ScenarioFileInputFormat,const char*,PetscInt*);
PETSC_EXTERN PetscErrorCode SOPFLOWGetObjective(SOPFLOW,PetscReal*);
PETSC_EXTERN PetscErrorCode SOPFLOWGetSolution(SOPFLOW,PetscInt,Vec*);
PETSC_EXTERN PetscErrorCode SOPFLOWPrintSolution(SOPFLOW,PetscInt);
PETSC_EXTERN PetscErrorCode SOPFLOWSaveSolution(SOPFLOW,PetscInt,OutputFormat,const char[]);
PETSC_EXTERN PetscErrorCode SOPFLOWSaveSolutionAll(SOPFLOW,OutputFormat,const char[]);
PETSC_EXTERN PetscErrorCode SOPFLOWGetConvergenceStatus(SOPFLOW,PetscBool*);
PETSC_EXTERN PetscErrorCode SOPFLOWGetMode(SOPFLOW,PetscInt*);
PETSC_EXTERN PetscErrorCode SOPFLOWGetNumIterations(SOPFLOW,PetscInt*);
PETSC_EXTERN PetscErrorCode SOPFLOWSetTolerance(SOPFLOW,PetscReal);
PETSC_EXTERN PetscErrorCode SOPFLOWGetTolerance(SOPFLOW,PetscReal*);

#endif


