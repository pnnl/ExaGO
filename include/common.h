/*
 * Public header file containing the common objects, API used by different applications
 */

#ifndef COMMON_H
#define COMMON_H

#include <petsc.h>

/** 
 * The communicator context 
 */
struct _p_COMM {
  MPI_Comm     type; /**< MPI communicator SELF or WORLD */
  PetscMPIInt  rank; /**< Process rank */
  PetscMPIInt  size; /**< Communicator size */
  PetscInt     refct; /**< Reference count to know how many objects are sharing the communicator */
};

typedef struct _p_COMM *COMM;

PETSC_EXTERN PetscErrorCode COMMCreate(MPI_Comm,COMM*);
PETSC_EXTERN PetscErrorCode COMMDestroy(COMM*);

#endif
