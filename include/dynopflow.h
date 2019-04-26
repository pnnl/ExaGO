/**
 * @file dynopflow.h
 * @brief Public header file containing the public API for dynamics optimization application
 *
 * Note: The functions on below are not defined on header file, but implemented on cpp file:\n
 * eval_dynopflow_f(PetscInt n, PetscScalar* x, Bool new_x, PetscScalar* obj_value, UserDataPtr user_data)\n
 * eval_dynopflow_grad_f(PetscInt n, PetscScalar* x, Bool new_x, PetscScalar* grad_f, UserDataPtr user_data)\n
 * eval_dynopflow_g(PetscInt n, PetscScalar* x, Bool new_x, PetscInt m, PetscScalar* g, UserDataPtr user_data)\n
 * eval_dynopflow_jac_g(PetscInt n, PetscScalar *x, Bool new_x, PetscInt m, PetscInt nele_jac, PetscInt *iRow, PetscInt *jCol, PetscScalar *values, UserDataPtr user_data)\n
 * eval_dynopflow_h(Index n, Number *x, Bool new_x, Number obj_factor, Index m, Number *lambda, Bool new_lambda, Index nele_hess, Index *iRow, Index *jCol, Number *values, UserDataPtr user_data)\n
 * DYNOPFLOWUpdateDYNParameters(DYNOPFLOW,Vec)\n
 * DYNOPFLOWComputeSensiP(DYNOPFLOW)\n
 * DYNOPFLOWComputeDYNJacobianP(TS ts,PetscReal t,Vec X,Mat jacP,void *ctx)\n
 * DYNOPFLOWCostIntegrand(TS ts,PetscReal t,Vec X,Vec R,void *ctx)\n
 * DYNOPFLOWDRDYFunction(TS ts,PetscReal t,Vec X,Vec *drdy,void *ctx)\n
 * DYNOPFLOWDRDPFunction(TS ts,PetscReal t,Vec U,Vec *drdp,void *ctx)\n
 * DYNOPFLOWCreateGlobalVector(DYNOPFLOW dynopflow,Vec *vec)\n
 * DYNOPFLOWSetVariableandConstraintBounds(DYNOPFLOW dynopflow,Vec Xl,Vec Xu,Vec Gl,Vec Gu)\n
 * DYNOPFLOWGetConstraintJacobianNonzeros(DYNOPFLOW dynopflow,PetscInt *nnz)\n
 * DYNOPFLOWComputeSensiP(DYNOPFLOW dynopflow)\n
 * DYNOPFLOWSetUp(DYNOPFLOW dynopflow)\n
 * DYNOPFLOWSolveInitialize(DYNOPFLOW dynopflow)\n
 */
#ifndef DYNOPFLOW_H
#define DYNOPFLOW_H

#include <ps.h>

typedef struct _p_DYNOPFLOW *DYNOPFLOW;

/**
 * @brief Creates a DYNOPFLOW application context
 * @param [in] MPI_Comm mpicomm - MPI communicator
 * @param [out] DYNOPFLOW* dynopflowout - the new DYNOPFLOW application context
 */
extern PetscErrorCode DYNOPFLOWCreate(MPI_Comm,DYNOPFLOW*);
/**
 * @brief Destroys the dynopflow application object that was created with DYNOPFLOWCreate().
 * @param [in] DYNOPFLOW* dynopflow - the DYNOPFLOW application context pointer
 */
extern PetscErrorCode DYNOPFLOWDestroy(DYNOPFLOW*);
/**
 * @brief Reads the network data given in MATPOWER data format 
 * @param [in] DYNOPFLOW DYNOPFLOW - The DYNOPFLOW object
 * @param [in] const char[] netfile - The name of the network file
 * Notes: With multiple processes, each process reads the data file.
 */
extern PetscErrorCode DYNOPFLOWReadMatPowerData(DYNOPFLOW,const char[]);
/**
 * @brief Reads the machine dynamic data file
 * @param [in] DYNOPFLOW DYNOPFLOW - The DYNOPFLOW object
 * @param [in] const char[] netfile - The name of the network file
 * Notes: With multiple processes, each process reads the data file.
 */
extern PetscErrorCode DYNOPFLOWReadDyrData(DYNOPFLOW,const char[]);
/**
 * @brief Reads the machine dynamic data file
 * @param [in] DYNOPFLOW DYNOPFLOW - The DYNOPFLOW object
 * @param [in] const char[] eventfile - the name of the event file
 * Notes: With multiple processes, each process reads the data file.
 */
extern PetscErrorCode DYNOPFLOWReadEventData(DYNOPFLOW,const char[]);
/**
 * @brief Sets the duration of the dynamic simulation
 * @param [in] DYNOPFLOW dynopflow - the dynopf application object
 * @param [in] PetscInt max_steps  - the maximum number of steps that the time-stepping solver takes
 * @param [in] PetscReal tend - the end time
 */
extern PetscErrorCode DYNOPFLOWSetDuration(DYNOPFLOW,PetscInt,PetscReal);
/**
 * @brief Sets the start time and the time step of the dynamic simulation
 * @param [in] DYNOPFLOW dynopflow - the dynopf application object
 * @param [in] PetscInt start_time - start time
 * @param [in] PetscReal time_step  - the step size (in seconds)
 * Notes: For variable time-stepping methods, this step is used as the initial time step.
 */
extern PetscErrorCode DYNOPFLOWSetStartTimeAndStep(DYNOPFLOW,PetscReal,PetscReal);
/**
 * @brief Solves the dynamics constrained optimal power flow
 * @param [in] DYNOPFLOW dynopflow - the dynopf application object
 */
extern PetscErrorCode DYNOPFLOWSolve(DYNOPFLOW);
#endif
