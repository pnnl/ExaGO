/**
 * @file common.h
 * @brief Public header file containing the common objects, API used by different applications
 */

#ifndef COMMON_H
#define COMMON_H

#include <petsc.h>

/** 
 * @brief The communicator context 
 */
struct _p_COMM {
  MPI_Comm     type; /**< MPI communicator SELF or WORLD */
  PetscMPIInt  rank; /**< Process rank */
  PetscMPIInt  size; /**< Communicator size */
  PetscInt     refct; /**< Reference count to know how many objects are sharing the communicator */
};

typedef struct _p_COMM *COMM;

/**
 * @brief Creates the communicator object COMM
 * @param [in] MPI_Comm mpicomm - The MPI communicator
 * @param [out] COMM* outcomm - The COMM object
 */
extern PetscErrorCode COMMCreate(MPI_Comm,COMM*);
/**
 * @brief Destroys the communicator object COMM created with COMMCreate
 * @param [in] COMM* outcomm - The COMM object
 */
extern PetscErrorCode COMMDestroy(COMM*);
/**
 * @brief NOT IMPLIMENTED
 */
extern PetscErrorCode SetMatrixValues(Mat,PetscInt,PetscInt[],PetscInt,PetscInt[],PetscScalar[]);
#endif
