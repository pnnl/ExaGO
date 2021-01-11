#ifndef EXAGO_VERSION_H
#define EXAGO_VERSION_H

#include <petsc.h>

/**
 * Print information about the current ExaGO build. The user is responsible for
 * free'ing this memory.
 */
PETSC_EXTERN PetscErrorCode ExaGOVersionGetFullVersionInfo(char**);

/**
 * Print release date for the current ExaGO build in format "%Y-%m-%d".
 * The user is responsible for free'ing this memory.
 */
PETSC_EXTERN PetscErrorCode ExaGOVersionGetReleaseDate(char**);

/**
 * Sets major, minor, and patch versions for current ExaGO build. The user is
 * responsible for free'ing this memory.
 */
PETSC_EXTERN PetscErrorCode ExaGOVersionGetVersion(int*,int*,int*);

/**
 * Sets string with build version for current ExaGO build in format
 * "major.minor.patch". The user is responsible for free'ing this memory.
 */
PETSC_EXTERN PetscErrorCode ExaGOVersionVersionStr(char**);

/** Number of ExaGO dependencies that are explicitly tracked */
#define ExaGONumDependencies 6

/**
 * Contains statically determinable information about the current ExaGO build.
 */
PETSC_EXTERN const char* ExaGODependencyNames[ExaGONumDependencies];

/**
 * Each entry indicates whether the dependency indicated by the same index
 * of array `ExaGODependencyNames` is enabled or not.
 */
PETSC_EXTERN const PetscBool ExaGOIsDependencyEnabled[ExaGONumDependencies];

#endif
