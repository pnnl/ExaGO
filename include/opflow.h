/**
 * @file opflow.h
 * @brief Public header file for optimal power flow application.
 *
 * Note: The functions on below are not defined on header file, but implemented on cpp file:\n
 * eval_opflow_f(PetscInt n, PetscScalar* x, Bool new_x, PetscScalar* obj_value, UserDataPtr user_data)\n
 * eval_opflow_grad_f(PetscInt n, PetscScalar* x, Bool new_x, PetscScalar* grad_f, UserDataPtr user_data)\n
 * eval_opflow_g(PetscInt n, PetscScalar* x, Bool new_x, PetscInt m, PetscScalar* g, UserDataPtr user_data)\n
 * eval_opflow_jac_g(PetscInt n, PetscScalar *x, Bool new_x, PetscInt m, PetscInt nele_jac, PetscInt *iRow, PetscInt *jCol, PetscScalar *values, UserDataPtr user_data)\n
 * OPFLOWGetLagrangianHessianNonzeros(OPFLOW opflow,PetscInt *nnz)\n
 * OPFLOWSetLagrangianHessianLocations(OPFLOW opflow, PetscInt *row, PetscInt *col)\n
 * OPFLOWSetLagrangianHessianValues(OPFLOW opflow, PetscScalar obj_factor, Vec X, Vec Lambda, PetscScalar *values)\n
 * OPFLOWCreateConstraintJacobian(OPFLOW opflow,Mat *mat)\n
 * OPFLOWGetConstraintJacobianNonzeros(OPFLOW opflow,PetscInt *nnz)\n
 * OPFLOWObjectiveFunction(OPFLOW opflow,Vec X, PetscScalar* obj)\n
 * OPFLOWObjGradientFunction(OPFLOW opflow,Vec X, Vec grad)\n
 * OPFLOWConstraintFunction(OPFLOW opflow,Vec X,Vec G)\n
 * OPFLOWSetConstraintJacobianLocations(OPFLOW opflow, PetscInt *row, PetscInt *col)\n
 * OPFLOWSetConstraintJacobianValues(OPFLOW opflow, Vec X,PetscScalar *values)\n
 */

#ifndef OPFLOW_H
#define OPFLOW_H

#include <ps.h>

typedef struct _p_OPFLOW *OPFLOW;

/**
 * @brief Creates an optimal power flow application object
 * @param [in] MPI_Comm mpicomm - The MPI communicator
 * @param [out] OPFLOW* opflowout - The optimal power flow application object
 */
extern PetscErrorCode OPFLOWCreate(MPI_Comm,OPFLOW*);
/**
 * @brief Destroys the optimal power flow application object
 * @param [in] OPFLOW* opflowout - The optimal power flow application object
 */
extern PetscErrorCode OPFLOWDestroy(OPFLOW*);
/**
 * @brief Reads the network data given in MATPOWER data format 
 * @param [in] OPFLOW opflow - The OPFLOW object
 * @param [in] const char[] netfile - The name of the network file
 */
extern PetscErrorCode OPFLOWReadMatPowerData(OPFLOW,const char[]);
/**
 * @brief Sets up a power flow application object
 * @param [in] OPFLOW opflowout - The optimal power flow application object
 * Notes:
 * This routine sets up the OPFLOW object and the underlying PS object. It
 * also distributes the PS object when used in parallel.
 */
extern PetscErrorCode OPFLOWSetUp(OPFLOW);
/**
 * @brief Returns a global vector of the appropriate size and distribution conforming to the distribution of the PS object.
 * @param [in] OPFLOW opflowout - The optimal power flow application object
 * @param [out] Vec* vec - the global vector
 * Notes:
 * OPFLOWSetUp() must be called before calling this routine.
 */
extern PetscErrorCode OPFLOWCreateGlobalVector(OPFLOW,Vec*);
/**
 * @brief NOT IMPLIMENTED
 */
extern PetscErrorCode OPFLOWCreateMatrix(OPFLOW,Mat*);
/**
 * @brief Solves the AC optimal power flow
 * @param [in] OPFLOW opflow - The OPFLOW object
 */
extern PetscErrorCode OPFLOWSolve(OPFLOW);
#endif


