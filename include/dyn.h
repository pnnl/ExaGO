/**
 * @file dyn.h
 * @brief Public header file for dynamic simulation application
 *
 * Note: The functions on below are not defined on header file, but implemented on cpp file:\n
 * DYNGetTS(DYN dyn, TS ts)\n
 * DYNGetPS(DYN dyn, PS ps)\n
 * DYNViewer(TS ts)\n
 * DYNSetPostStepCallback(DYN dyn,PetscErrorCode (*poststepcallback)(DYN dyn))\n
 * DYNUserMonitor(TS ts,PetscInt stepnum,PetscReal t,Vec X,void *ctx)\n
 * CheckAlgebraicPart(DYN dyn,Vec v)\n
 * DYNEventMonitor_Machine(DYNEvent dynevent, PetscReal t, Vec X, PetscScalar *fval)\n
 * DYNEventPostFunction_Machine(DYNEvent dynevent, PetscInt nmonitors, PetscInt monidx[], PetscReal t, Vec localX, PetscBool forwardsolve, PetscBool *solve_alg)\n
 * DYNEventPostDGamma_Machine(DYNEvent dynevent, PetscInt nmonitors, PetscInt monidx[], PetscReal t, Vec localX, PetscBool forwardsolve)\n
 * DYNCheckandSetMachineLimits(DYN dyn, PetscReal t, Vec X,PetscBool forwardsolve)\n
 * DYNSetMachineEvents(DYN dyn)\n
 * DYNAdjointMonitor(TS ts,PetscInt stepnum,PetscReal t,Vec X,PetscInt ncost,Vec *lambda,Vec *mu,void *ctx)\n
 * DYNSetGenDispatchandStatus(DYN dyn, PetscInt busnum, PetscInt gennum, PetscInt status, PetscScalar Pg, PetscScalar Qg)\n
 * DYNComputePrePflow(DYN dyn)\n
 * DYNComputeSensiP(DYN dyn)\n
 * DYNComputeAdjointJacobianP(TS ts,PetscReal t,Vec X,Mat jacP,void *ctx)\n
 * DYNCostIntegrand(TS ts,PetscReal t,Vec X,Vec R,void *ctx)\n
 * DYNDRDYFunction(TS ts,PetscReal t,Vec X,Vec *drdy,void *ctx)\n
 * DYNDRDPFunction(TS ts,PetscReal t,Vec U,Vec *drdp,void *ctx)\n
 * DYNIJacobian(TS ts,PetscReal t,Vec X,Vec Xdot,PetscReal a, Mat J, Mat Jpre, void *ctx)\n
 * DYNIFunction(TS ts,PetscReal t,Vec X,Vec Xdot,Vec F,void *ctx)\n
 * DYNSetInitialConditions(DYN dyn,Vec X)\n
 * DYNSetLambdaToUVectors(DYN dyn,Vec* lambda,PetscInt ncostfcns)\n
 * DYNGetCostFunction(DYN dyn,PetscReal t,PetscReal *costfcn,PetscInt ncostfun) //Comments in .c file is also wrong\n
 * DYNSetUpDifferentialEqIS(DYN dyn, IS *isdiff)\n
 * DYNComputeForwardJacobianP(TS ts,PetscReal t,Vec X,Vec *jacP,void *ctx)\n
 */

#ifndef DYN_H
#define DYN_H

#include <ps.h>
#include <pflow.h>

typedef struct _p_DYN *DYN; /**< The dynamics simulation object */

/**
 * @brief Creates a dynamic simulation application object
 * @param [in] MPI_Comm mpicomm - The MPI communicator
 * @param [out] DYN* DYN- The dynamic simulation application object
 */
PETSC_EXTERN PetscErrorCode DYNCreate(MPI_Comm,DYN*);
/**
 * @brief Destroys the DYN object
 * @param [in] DYN* dyn- The dynamic simulation application object
 */
PETSC_EXTERN PetscErrorCode DYNDestroy(DYN*);
/**
 * @brief Sets up a DYN object
 * @param [in] DYN dyn- The dynamic simulation application object
 * Notes:
 * This routine sets up the DYN object and the underlying PS object. It
 * also distributes the PS object when used in parallel.
 */
PETSC_EXTERN PetscErrorCode DYNSetUp(DYN);
/**
 * @brief Runs the dynamic simulation
 * @param [in] DYN dyn- The dynamic simulation application object
 */
PETSC_EXTERN PetscErrorCode DYNSolve(DYN);
/**
 * @brief Reads the network data given in MATPOWER data format
 * @param [in] DYN DYN - The DYN object
 * @param [in] const char[] netfile - The name of the network file
 */
PETSC_EXTERN PetscErrorCode DYNReadMatPowerData(DYN, const char[]);
/**
 * @brief Reads the network data given in PSSE raw data format
 * @param [in] DYN DYN - The DYN object
 * @param [in] const char[] netfile - The name of the network file
 */
PETSC_EXTERN PetscErrorCode DYNReadPSSERawData(DYN, const char[]);
/**
 * @brief Reads the data file with dynamic models
 * @param [in] DYN DYN - The DYN object
 * @param [in] const char[] dyrfile - The name of the dyr file
 */
PETSC_EXTERN PetscErrorCode DYNReadDyrData(DYN, const char[]);
/**
 * @brief NOT IMPLIMENTED
 */
PETSC_EXTERN PetscErrorCode DYNReadEventData(DYN,const char[]);
/**
 * @brief Sets the duration of the dynamic simulation
 * @param [in] DYN DYN - The DYN object
 * @param [in] PetscInt max_steps  - the maximum number of steps that the time-stepping solver takes
 * @param [in] PetscReal tend - the end time
 */
PETSC_EXTERN PetscErrorCode DYNSetDuration(DYN, PetscInt,PetscReal);
/**
 * @brief Sets the start time and the time step of the dynamic simulation
 * @param [in] DYN DYN - The DYN object
 * @param [in] PetscReal start_time - start time
 * @param [in] PetscReal time_step  - the step size (in seconds)
 * Notes: For variable time-stepping methods, this step is used as the initial time step.
 */
PETSC_EXTERN PetscErrorCode DYNSetStartTimeAndStep(DYN,PetscReal,PetscReal);
/**
 * @brief Returns a global vector of the appropriate size and distribution conforming to the distribution of the PS object.
 * @param [in] DYN DYN - The DYN object
 * @param [out] Vec* vec - the global vector
 * Notes: DYNSetUp() must be called before calling this routine.
 */
PETSC_EXTERN PetscErrorCode DYNCreateGlobalVector(DYN, Vec*);
/**
 * @brief Returns a distributed matrix of appropriate size that can be used as the Jacobian
 * @param [in] DYN DYN - The DYN object
 * @param [out] Mat* mat - the matrix
 * Notes: DYNSetUp() must be called before calling this routine.
 */
PETSC_EXTERN PetscErrorCode DYNCreateMatrix(DYN,Mat*);
/**
 * @brief Sets up a DYN object
 * @param [in] DYN dyn- The dynamic simulation application object
 * Notes:
 * This routine sets up the DYN object and the underlying PS object. It
 * also distributes the PS object when used in parallel
 */
PETSC_EXTERN PetscErrorCode DYNAdjointSetUp(DYN);
/**
 * @brief Runs the dynamic simulation
 * @param [in] DYN dyn- The dynamic simulation application object
 * Notes:
 * This routine sets up the DYN object and the underlying PS object. It
 * also distributes the PS object when used in parallel.
 */
PETSC_EXTERN PetscErrorCode DYNAdjointSolve(DYN);
/**
 * @brief Destroys the DYN object
 * @param [in] DYN* DYN- The dynamic simulation application object
 */
PETSC_EXTERN PetscErrorCode DYNAdjointDestroy(DYN*);
#ifdef FWDSA
/**
 * @brief Sets up a DYN object
 * @param [in] DYN* DYN- The dynamic simulation application object
 * Notes: This routine sets up the DYN object and the underlying PS object. It
 * also distributes the PS object when used in parallel.
 */
PETSC_EXTERN PetscErrorCode DYNForwardSetUp(DYN);
/**
 * @brief Perform forward sensitivity analysis of the dynamic simulation
 * @param [in] DYN DYN- The dynamic simulation application object
 * Notes: The parameters considered are the bus voltage magnitudes and angles, and the generator real and reactive power
 * output. For each bus, its parameters are ordered as [Va, Vm, PG1, QG1, PG2, QG2,.., PGn, QGn]
 */
PETSC_EXTERN PetscErrorCode DYNForwardSolve(DYN);
/**
 * @brief Destroys the DYN object
 * @param [in] DYN* DYN- The dynamic simulation application object
 */
PETSC_EXTERN PetscErrorCode DYNForwardDestroy(DYN*);
/**
 * @brief Solve the algebraic part to obtain consistent sensitivity variables
 * @param [Unknown] SNES tssnes
 * @param [in] DYN* DYN- The dynamic simulation application object
 * @param [Unknown] PetscRea t
 * @param [Unknown] Vec X
 */
PETSC_EXTERN PetscErrorCode DYNForwardAlgebraicSolve(SNES,DYN,PetscReal,Vec);
#endif
/**
 * @brief Initial power flow solution for DYN
 * @param [in] DYN dyn - the DYN object
 * @param [in] PetscBool flag - compute initial power flow (0 = No, 1 = Yes)
 * Notes:
 * Calling this routine with PETSC_TRUE flag computes the steady-state power flow solution for the DYN object.
 */
PETSC_EXTERN PetscErrorCode DYNSetComputePrePflow(DYN,PetscBool);
/**
 * @brief This function runs the forward solve and returns the cost functions along with their sensitivities
 * @param [in] int argc
 * @param [in] char** **argv
 * @param [in] const char[] netfile - Name of the network file
 * @param [in] const char[] dyrfile - Name of the dynamic data dyr file
 * @param [in] const char[] eventfile - Name of the event file
 * @param [in] double* costfcnout - A double array having the cost functions (must be allocated)
 * @param [in] double* sensitivitiesout - An array holding the sensitivities (must be allocated)
 * @param [in] PetscBool single_costfcn - Flag to indicate whether a single
 * @param [in] PetscBool active_power_only - Flag to indicate that sensitivities w.r.t. only active power only needs to be computed.
 * @param [in] PetscBool updatedispatch - Flag to update dispatch
 * @param [in] PetscInt* genbus - An array of generator bus numbers
 * @param [in] PetscInt* gennum - An array of generator numbers (starting from 0)
 * @param [in] PetscInt* gen_status - An array of generator commitments
 * @param [in] PetscScalar* Pg - An array of active power dispatch
 * @param [in] PetscScalar* Qg - An array of reactvie power dispatch
 * @param [in] PetscReal endtime - End time of the simulation
 * @param [in] PetscReal stepsize - Time step size
 * @param [in] PetscInt cost_type - Flag to print voltage
 * @param [in] PetscInt print_vm - Indicator of cost function type (frequency voilation or frequency)
 * @param [in] PetscReal - Minimum frequency threshold
 * @param [in] PetscReal - Maximum frequency threshold
 * @param [in] PetscBool - Indicator of whether to use a monitor
 * @param [in] PetscReal* - Final time at which the simulation actually ends
 */
PETSC_EXTERN PetscInt DYNGetCostFunctionAndSensitivities(int,char**,const char[],const char[],const char[],double*,double*,PetscBool,PetscBool,PetscBool,PetscInt*,PetscInt*,PetscInt*,PetscScalar*,PetscScalar*,PetscReal,PetscReal,PetscInt,PetscInt,PetscReal,PetscReal,PetscBool,PetscReal*);
/**
 * @brief NO DESCRIPTION
 */
PETSC_EXTERN PetscErrorCode  Initialize(int ,char **);
/**
 * @brief NO DESCRIPTION
 */
PETSC_EXTERN PetscErrorCode  Finalize();
/**
 * @brief NOT IMPLIMENTED
 */
PETSC_EXTERN PetscErrorCode DYNSetGenDispatchandStatus(DYN,PetscInt,PetscInt,PetscInt,PetscScalar,PetscScalar);

/**
 * @brief Set Generator status
 * @param [in] DYN - The dynamics simulation object
 * @param [in] gbus - Generator bus number
 * @param [in] gid - Generator id
 * @param [in] status - Generator status (1 = ON, 0 = OFF)
 */
PETSC_EXTERN PetscErrorCode DYNSetGenStatus(DYN,PetscInt,const char*,PetscInt);
/**
 * @brief Perform finite difference of the dynamic simulation
 * @param [in] DYN DYN- The dynamic simulation application object
 * Notes: The parameters considered are the bus voltage magnitudes and angles, and the generator real and reactive power
 * output. For each bus, its parameters are ordered as [Va, Vm, PG1, QG1, PG2, QG2,.., PGn, QGn]
 */
PETSC_EXTERN PetscErrorCode DYNFiniteDiffSolve(DYN);

PETSC_EXTERN PetscErrorCode DYNSetPreStepCallback(DYN,PetscErrorCode(*)(DYN));
PETSC_EXTERN PetscErrorCode DYNSetPostStepCallback(DYN,PetscErrorCode (*)(DYN));
/**
 * @brief Sets a callback function for every successful step taken by DYN
 * @param [in] DYN - The dynamic simulation object
 * @param [in] poststepcallback - Callback function called after each successful time-step
 */

PETSC_EXTERN PetscErrorCode DYNGetTS(DYN,TS*);
/**
 * @brief Returns the PETSc time-stepping solver object used with DYN
 * @param [in] DYN - The dynamic simulation object
 * @param [out] ts - the PETSc time-stepping solver object
 */

PETSC_EXTERN PetscErrorCode DYNGetPS(DYN,PS*);
/**
 * @brief Returns the power system object PS used with DYN
 * @param [in] DYN - The dynamic simulation object
 * @param [out] ps - the power system object
 */

/**
 * @brief Returns the initial power flow (used to calculate PF at t=0) used with DYN
 * @param [in] DYN - The dynamic simulation object
 * @param [out] pflow - power flow object used in calculating initial conditions.
 */

PETSC_EXTERN PetscErrorCode DYNGetInitPflow(DYN,PFLOW*);
PETSC_EXTERN PetscErrorCode DYNComputePrePflow(DYN,PetscBool*);
PETSC_EXTERN PetscErrorCode DYNSetUserData(DYN,void*);
PETSC_EXTERN PetscErrorCode DYNGetUserData(DYN,void**);
PETSC_EXTERN PetscErrorCode DYNGetTime(DYN,PetscReal*);
PETSC_EXTERN PetscErrorCode DYNGetStepNumber(DYN,PetscInt*);
PETSC_EXTERN PetscErrorCode DYNSaveSolution(DYN,PetscReal,Vec);
PETSC_EXTERN PetscErrorCode DYNGetBusVoltage(DYN,PetscInt,PetscScalar*,PetscScalar*,PetscBool*);
PETSC_EXTERN PetscErrorCode DYNSetBusVoltage(DYN,PetscInt,PetscScalar,PetscScalar);
PETSC_EXTERN PetscErrorCode DYNGetBusGenCurrent(DYN,PetscInt,PetscScalar*,PetscScalar*,PetscBool*);
PETSC_EXTERN PetscErrorCode DYNSetBusConstantPowerLoad(DYN,PetscInt,PetscScalar,PetscScalar);
PETSC_EXTERN PetscErrorCode DYNSolveAlgebraicEquationsOnly(DYN);
PETSC_EXTERN PetscErrorCode DYNIsSwitchingSolutionStep(DYN,PetscBool*);
PETSC_EXTERN PetscErrorCode DYNInitializeEvents(DYN);

PETSC_EXTERN PetscErrorCode DYNCheckandSetMachineLimits(DYN,PetscReal,Vec,PetscBool);

PETSC_EXTERN PetscErrorCode DYNGetCurrentStepSolution(DYN,Vec*);

PETSC_EXTERN PetscErrorCode DYNSetTimeStep(DYN,PetscReal);
PETSC_EXTERN PetscErrorCode DYNGetTimeStep(DYN,PetscReal*);
PETSC_EXTERN PetscErrorCode DYNRollbackStep(DYN);

#endif
