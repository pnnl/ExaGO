#ifndef EXAGO_UTILS_H
#define EXAGO_UTILS_H
#include <petsc.h>

typedef enum {
  EXAGO_LOG_INFO=0,
  EXAGO_LOG_WARN,
  EXAGO_LOG_ERROR,
  EXAGO_LOG_DISABLE,
  EXAGO_LOG_NUM_LOG_LEVELS,
} ExaGOVerbosityLevel;

PETSC_EXTERN const char* ExaGOVerbosityNames[EXAGO_LOG_NUM_LOG_LEVELS];

/**
 * Set the name for the logfile to be used; must be set before
 * `ExaGOInitialize` is called.
 */
PETSC_EXTERN PetscErrorCode ExaGOLogSetLoggingFileName(char*);

/** Retrieve the filename being used for the logfile */
PETSC_EXTERN PetscErrorCode ExaGOLogGetLoggingFileName(char**);

/** Indicates whether a logfile is being used or stdout/stderr */
PETSC_EXTERN PetscErrorCode ExaGOLogIsUsingLogFile(PetscBool*);

/** Get minimum loglevel for a log to be printed */
PETSC_EXTERN PetscErrorCode ExaGOLogGetMinLogLevel(ExaGOVerbosityLevel*);

/** Set minimum loglevel for a log to be printed */
PETSC_EXTERN PetscErrorCode ExaGOLogSetMinLogLevel(ExaGOVerbosityLevel);

#if !defined(EXAGO_DISABLE_LOGGING)

/** 
 * @brief Implementation to log string according to ExaGO build configuration.
 * @note Users should use `ExaGOLog` instead, as this can be disabled for
 * optimization.
 */
PETSC_EXTERN void ExaGOLogImpl(ExaGOVerbosityLevel,char*,...);

/** `ExaGOLog` is the user-facing logging function */
#define ExaGOLog(level,fmt,...) ExaGOLogImpl(level,fmt,__VA_ARGS__)

#else

/** If logging is disabled, logging is a no-op */
#define ExaGOLog(x,y,...) (void)(x)

#endif

PETSC_EXTERN PetscErrorCode ExaGOLogIsInitialized(PetscBool*);

/**
 * Initialize an ExaGO application.
 *
 * @note this takes care of Petsc initialization, so don't this function in
 * conjunction with `PetscInitialize.`
 */
PETSC_EXTERN PetscErrorCode ExaGOInitialize(MPI_Comm,int*argc,char***argv,
    char*appname,char*help);

/**
 * Teardown for an ExaGO application.
 *
 * @note this takes care of Petsc finalization, so don't this function in
 * conjunction with `PetscFinalize`.
 */
PETSC_EXTERN PetscErrorCode ExaGOFinalize();

/** Returns 1 if files exists, else 0 */
PETSC_EXTERN int doesFileExist(char*);

/** Returns 1 if directory exists, else 0 */
PETSC_EXTERN int doesDirExist(char*);

/** first i in [0, npths) such that pths[i] is statable as a regular file */
PETSC_EXTERN int anyFileExist(char**, int);

/** Determines if scaled difference between two reals is within tolerance */
PETSC_EXTERN int isEqual(double,double,double tol,double*err);

/** Determines if difference between two ints is within tolerance */
PETSC_EXTERN int isEqualInt(int value,int reference,int tol,int*err);

#endif
