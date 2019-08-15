/**
 * @file scopflow.h
 * @brief Public header file for security constrained optimal power flow application.
 *
 */

#ifndef SCOPFLOW_H
#define SCOPFLOW_H

#include <opflow.h>

typedef struct _p_SCOPFLOW *SCOPFLOW;

/**
 * @brief Creates an security constrained optimal power flow application object
 * @param [in] MPI_Comm mpicomm - The MPI communicator
 * @param [out] OPFLOW* opflowout - The scopf application object
 */
PETSC_EXTERN PetscErrorCode SCOPFLOWCreate(MPI_Comm,SCOPFLOW*);
/**
 * @brief Destroys the scopf application object
 * @param [in] SCOPFLOW* scopflow - The scopflow application object
 */
PETSC_EXTERN PetscErrorCode SCOPFLOWDestroy(SCOPFLOW*);
/**
 * @brief Reads the network data given in MATPOWER data format 
 * @param [in] SCOPFLOW scopflow - The SCOPFLOW object
 * @param [in] const char[] netfile - The name of the network file
 */
PETSC_EXTERN PetscErrorCode SCOPFLOWSetNetworkData(SCOPFLOW,const char[]);
/**
 * @brief Sets up a security constrained power flow application object
 * @param [in] SCOPFLOW scopflow - The security constrained optimal power flow application object
 * Notes:
 * This routine sets up the SCOPFLOW object and the underlying PS object.
 */
PETSC_EXTERN PetscErrorCode SCOPFLOWSetUp(SCOPFLOW);

/**
 * @brief Returns a global vector of the appropriate size and distribution conforming to the distribution of the PS object.
 * @param [in] SCOPFLOW scopflow - The security constrained optimal power flow application object
 * @param [out] Vec* vec - the global vector
 * Notes:
 * SCOPFLOWSetUp() must be called before calling this routine.
 */
PETSC_EXTERN PetscErrorCode SCOPFLOWCreateGlobalVector(SCOPFLOW,Vec*);
/**
 * @brief NOT IMPLIMENTED
 */
PETSC_EXTERN PetscErrorCode SCOPFLOWCreateMatrix(SCOPFLOW,Mat*);
/**
 * @brief Solves the AC SCOPF optimal power flow
 * @param [in] OPFLOW opflow - The OPFLOW object
 */
PETSC_EXTERN PetscErrorCode SCOPFLOWSolve(SCOPFLOW);

PETSC_EXTERN PetscErrorCode SCOPFLOWSetNumScenarios(SCOPFLOW,PetscInt);

#endif


