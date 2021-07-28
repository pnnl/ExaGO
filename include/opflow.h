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
#define OPFLOWSOLVER_HIOPSPARSE "HIOPSPARSE"

typedef struct _p_OPFLOW *OPFLOW;

typedef enum {
  MIN_GEN_COST, /* Generator cost minimization */
  MIN_GENSETPOINT_DEVIATION, /* Minimize generator set-point deviation */
  NO_OBJ, /* No objective function */
}OPFLOWObjectiveType;

/* Generator bus voltage type */
typedef enum {
  VARIABLE_WITHIN_BOUNDS, /* Voltage can vary with voltage bounds */
  FIXED_WITHIN_QBOUNDS,  /* Voltage is fixed within reactive power bounds */
  FIXED_AT_SETPOINT, /* Voltage is fixed to the set-point voltage */
}OPFLOWGenBusVoltageType;

typedef enum {
  OPFLOWINIT_MIDPOINT, /* Midpoint */
  OPFLOWINIT_FROMFILE, /* From file */
  OPFLOWINIT_ACPF,      /* From AC power flow solution */
  OPFLOWINIT_FLATSTART  /* Voltage flat start */
}OPFLOWInitializationType; 

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
PETSC_EXTERN PetscErrorCode OPFLOWSetTolerance(OPFLOW,PetscReal);
PETSC_EXTERN PetscErrorCode OPFLOWGetTolerance(OPFLOW,PetscReal*);
PETSC_EXTERN PetscErrorCode OPFLOWGetNumIterations(OPFLOW,PetscInt*);
PETSC_EXTERN PetscErrorCode OPFLOWGetObjectiveType(OPFLOW,OPFLOWObjectiveType*);
PETSC_EXTERN PetscErrorCode OPFLOWSetObjectiveType(OPFLOW,OPFLOWObjectiveType);
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

PETSC_EXTERN PetscErrorCode OPFLOWGetVariableOrdering(OPFLOW,int**);
PETSC_EXTERN PetscErrorCode OPFLOWGetSizes(OPFLOW,int*,int*,int*);

PETSC_EXTERN PetscErrorCode OPFLOWHasGenSetPoint(OPFLOW,PetscBool);
PETSC_EXTERN PetscErrorCode OPFLOWHasLoadLoss(OPFLOW,PetscBool);
PETSC_EXTERN PetscErrorCode OPFLOWHasBusPowerImbalance(OPFLOW,PetscBool);
PETSC_EXTERN PetscErrorCode OPFLOWUseAGC(OPFLOW,PetscBool);

PETSC_EXTERN PetscErrorCode OPFLOWSetGenBusVoltageType(OPFLOW,OPFLOWGenBusVoltageType);

PETSC_EXTERN PetscErrorCode OPFLOWGetObjectiveType(OPFLOW,OPFLOWObjectiveType*);
PETSC_EXTERN PetscErrorCode OPFLOWSetObjectiveType(OPFLOW,OPFLOWObjectiveType);

PETSC_EXTERN PetscErrorCode OPFLOWSetInitializationType(OPFLOW, OPFLOWInitializationType);

/* OPFLOWGetPS - Gets the underlying PS object

  Input Parameters                                                            
. opflow - the OPFLOW object                  

  Output Parameters             
. ps - the ps object       

  Notes: This function returns the PS object that holds the network data. Using the PS object            
         one can make changes to the network parameters. A typical case is
         changing some network parameters before solving opflow.                      
         OPFLOWSetUpPS() must be called before OPFLOWGetPS()
*/
PETSC_EXTERN PetscErrorCode OPFLOWGetPS(OPFLOW,PS*);

/* OPFLOWSetUpPS - Sets the underlying PS network object to be used by OPFLOW

  Input Parameters                                                             
. opflow - the OPFLOW object                                                                        
 
  Notes: This function is an intermediate function that can be called for setting up
         the PS network object prior to solving OPFLOW. A typical use-case is some
         network parameter needs changing before solving opflow. In such case,
         the work flow would be                                                                                
  1. OPFLOWCreate();                                                                                                 
  2. OPFLOWReadMatPowerData();
  3. OPFLOWSetUpPS();           
  4. OPFLOWGetPS();             
  ... change the network data by the PS object retrieved via OPFLOWGetPS().
  5. OPFLOWSolve();                                                                                   
 Skip steps 3 and 4 if no network changes are needed to be done.             
*/
PETSC_EXTERN PetscErrorCode OPFLOWSetUpPS(OPFLOW);

/*
  OPFLOWSolutionToPS - Updates the PS struct from OPFLOW solution

  Input Parameters:
. opflow - the OPFLOW object

  Notes: Updates the different fields in the PS struct from the OPFLOW solution.
         This function should be called after OPFLOWSolve() before retrieving
	 values from the PS struct.
*/
PETSC_EXTERN PetscErrorCode OPFLOWSolutionToPS(OPFLOW);

typedef PetscErrorCode (*OPFLOWAuxObjectiveFunction)(OPFLOW,const double*,double*,void*);
typedef PetscErrorCode (*OPFLOWAuxGradientFunction)(OPFLOW,const double*,double*,void*);
typedef PetscErrorCode (*OPFLOWAuxHessianFunction)(OPFLOW,const double*,Mat,void*);

PETSC_EXTERN PetscErrorCode OPFLOWSetAuxillaryObjective(OPFLOW,OPFLOWAuxObjectiveFunction,OPFLOWAuxGradientFunction,OPFLOWAuxHessianFunction,void*);

PETSC_EXTERN PetscErrorCode OPFLOWSetUpdateVariableBoundsFunction(OPFLOW,PetscErrorCode (*)(OPFLOW,Vec,Vec,void*),void*);

#endif


