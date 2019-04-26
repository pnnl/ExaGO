/**
 * @file pflow.h
 * @brief Public header file for power flow application
 *
 * Note: The functions on below are not defined on header file, but implemented on cpp file:\n
 * PFLOWSetInitialGuess(PFLOW pflow,Vec X)\n
 * PFLOWFunction(SNES snes,Vec X,Vec F,void *ctx)\n
 * PFLOWJacobian(SNES snes,Vec X,Mat J, Mat Jpre, void *ctx)\n
 */

#ifndef PFLOW_H
#define PFLOW_H

#include <ps.h>

typedef struct _p_PFLOW *PFLOW;

/**
 * @brief Creates a power flow application object
 * @param [in] MPI_Comm mpicomm - The MPI communicator
 * @param [out] PFLOW* pflowout - The power flow application object
 */
PETSC_EXTERN PetscErrorCode PFLOWCreate(MPI_Comm,PFLOW*);
/**
 * @brief Destroys the power flow application object
 * @param [in] PFLOW* pflow - The PFLOW object to destroy
 */
PETSC_EXTERN PetscErrorCode PFLOWDestroy(PFLOW*);
/**
 * @brief Reads the network data given in MATPOWER data format 
 * @param [in] PFLOW pflow - The PFLOW object to destroy
 * @param [in] const char[] netfile - The name of the network file
 */
PETSC_EXTERN PetscErrorCode PFLOWReadMatPowerData(PFLOW,const char[]);
/**
 * @brief Reads the network data given in PSSE raw data format 
 * @param [in] PFLOW pflow - The PFLOW object to destroy
 * @param [in] const char[] netfile - The name of the network file
 */
PETSC_EXTERN PetscErrorCode PFLOWReadPSSERawData(PFLOW,const char[]);
/**
 * @brief Sets up a power flow application object
 * @param [in] PFLOW pflow - The PFLOW object to destroy
 * Notes:
 * This routine sets up the PFLOW object and the underlying PS object. It
 * also distributes the PS object when used in parallel.
 */
PETSC_EXTERN PetscErrorCode PFLOWSetUp(PFLOW);
/**
 * @brief Returns a global vector of the appropriate size and distribution conforming to the distribution of the PS object.
 * @param [in] PFLOW pflow - The PFLOW object to destroy
 * @param [out] Vec* vec - the global vector
 * Notes:
 * PFLOWSetUp() must be called before calling this routine.
 */
PETSC_EXTERN PetscErrorCode PFLOWCreateGlobalVector(PFLOW,Vec*);
/**
 * @brief Returns a distributed matrix of appropriate size that can be used as the Jacobian
 * @param [in] PFLOW pflow - The PFLOW object to destroy
 * @param [out] Mat* mat - the matrix
 * Notes:
 * PFLOWSetUp() must be called before calling this routine.
 */
PETSC_EXTERN PetscErrorCode PFLOWCreateMatrix(PFLOW,Mat*);
/**
 * @brief Solves the AC power flow equations
 * @param [in] PFLOW pflow - The PFLOW object to destroy
 */
PETSC_EXTERN PetscErrorCode PFLOWSolve(PFLOW);
/**
 * @brief Updates the buses and the branches with the solution from the power flow
 * @param [in] PFLOW pflow - The PFLOW object to destroy
 */
PETSC_EXTERN PetscErrorCode PFLOWPostSolve(PFLOW);
PETSC_EXTERN PetscErrorCode PFLOWConverged(PFLOW,PetscBool*);

/**
 * @brief Sets the line status
 * @param [in] pflow - The PFLOW object
 * @param [in] fbus  - From bus
 * @param [in] tbus  - To bus
 * @param [in] id    - line id
 * @param [in] status - line status (0 = off, 1 = on)
 */
PETSC_EXTERN PetscErrorCode PFLOWSetLineStatus(PFLOW,PetscInt,PetscInt,const char*,PetscInt);
PETSC_EXTERN PetscErrorCode PFLOWSetGenStatus(PFLOW,PetscInt,const char*,PetscInt);

/**
 * @brief Returns the bus voltage magnitude and angle if bus found, otherwise sets the
 * found flag to FALSE.
 * @param [in] pflow - The PFLOW object
 * @param [in] bus  - the bus number
 * @param [in] Vm  -  voltage magnitude
 * @param [in] Va    - voltage angle (degrees)
 * @param [in] status - line status (0 = off, 1 = on)
 */
PETSC_EXTERN PetscErrorCode PFLOWGetBusVoltage(PFLOW,PetscInt,PetscScalar*,PetscScalar*,PetscBool*);
/**
 * @brief Sets the bus voltage magnitude and angle for a given bus
 * @param [in] pflow - The PFLOW object
 * @param [in] busnum  - the bus number
 * @param [in] Vm  -  voltage magnitude
 * @param [in] Va    - voltage angle (degrees)
 */
PETSC_EXTERN PetscErrorCode PFLOWSetBusVoltage(PFLOW,PetscInt,PetscScalar,PetscScalar);

PETSC_EXTERN PetscErrorCode PFLOWGetGenDispatch(PFLOW,PetscInt,PetscInt,PetscScalar*,PetscScalar*,PetscBool*);
PETSC_EXTERN PetscErrorCode PFLOWSetGenDispatch(PFLOW,PetscInt,PetscInt,PetscScalar,PetscScalar);

PETSC_EXTERN PetscErrorCode PFLOWSetLoadPower(PFLOW,PetscInt,PetscScalar,PetscScalar);
PETSC_EXTERN PetscErrorCode PFLOWGetLoadPower(PFLOW,PetscInt,PetscScalar*,PetscScalar*,PetscBool*);

PETSC_EXTERN PetscErrorCode PFLOWAddBusShunt(PFLOW,PetscInt,PetscScalar,PetscScalar);


#endif


