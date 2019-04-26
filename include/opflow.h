/**
 * @file opflow.h
 * @brief Public header file for optimal power flow application.
 *
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
PETSC_EXTERN PetscErrorCode OPFLOWCreate(MPI_Comm,OPFLOW*);
/**
 * @brief Destroys the optimal power flow application object
 * @param [in] OPFLOW* opflowout - The optimal power flow application object
 */
PETSC_EXTERN PetscErrorCode OPFLOWDestroy(OPFLOW*);
/**
 * @brief Reads the network data given in MATPOWER data format 
 * @param [in] OPFLOW opflow - The OPFLOW object
 * @param [in] const char[] netfile - The name of the network file
 */
PETSC_EXTERN PetscErrorCode OPFLOWReadMatPowerData(OPFLOW,const char[]);
/**
 * @brief Sets up a power flow application object
 * @param [in] OPFLOW opflowout - The optimal power flow application object
 * Notes:
 * This routine sets up the OPFLOW object and the underlying PS object. It
 * also distributes the PS object when used in parallel.
 */
PETSC_EXTERN PetscErrorCode OPFLOWSetUp(OPFLOW);
/**
 * @brief Returns a global vector of the appropriate size and distribution conforming to the distribution of the PS object.
 * @param [in] OPFLOW opflowout - The optimal power flow application object
 * @param [out] Vec* vec - the global vector
 * Notes:
 * OPFLOWSetUp() must be called before calling this routine.
 */
PETSC_EXTERN PetscErrorCode OPFLOWCreateGlobalVector(OPFLOW,Vec*);
/**
 * @brief NOT IMPLIMENTED
 */
PETSC_EXTERN PetscErrorCode OPFLOWCreateMatrix(OPFLOW,Mat*);
/**
 * @brief Solves the AC optimal power flow
 * @param [in] OPFLOW opflow - The OPFLOW object
 */
PETSC_EXTERN PetscErrorCode OPFLOWSolve(OPFLOW);
#endif


