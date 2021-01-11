#ifndef EXAGO_UTILSIMPL_H
#define EXAGO_UTILSIMPL_H
#include <utils.h>

/** Get file pointer for ExaGO logging */
PETSC_EXTERN PetscErrorCode ExaGOLogGetLoggingFilePointer(FILE**);

/** Set file pointer for ExaGO logging */
PETSC_EXTERN PetscErrorCode ExaGOLogSetLoggingFilePointer(FILE*);

/** Set communicator for ExaGO logging */
PETSC_EXTERN PetscErrorCode ExaGOLogSetComm(MPI_Comm);

/** Configures ExaGO logging system */
PETSC_EXTERN PetscErrorCode ExaGOLogInitialize();

/** Destroys ExaGO logging system */
PETSC_EXTERN PetscErrorCode ExaGOLogFinalize();

#endif
