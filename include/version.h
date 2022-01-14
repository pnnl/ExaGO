#ifndef EXAGO_VERSION_H
#define EXAGO_VERSION_H

#include <petsc.h>
#include <unordered_map>

/**
 * Print information about the current ExaGO build. The user is responsible for
 * free'ing this memory.
 */
PETSC_EXTERN PetscErrorCode ExaGOVersionGetFullVersionInfo(std::string &);

/**
 * Sets major, minor, and patch versions for current ExaGO build. The user is
 * responsible for free'ing this memory.
 */
PETSC_EXTERN PetscErrorCode ExaGOVersionGetVersion(int *, int *, int *);

/**
 * Sets string with build version for current ExaGO build in format
 * "major.minor.patch". The user is responsible for free'ing this memory.
 */
PETSC_EXTERN PetscErrorCode ExaGOVersionGetVersionStr(std::string &);

/**
 *  @brief Returns map with keys as the names of the dependencies and values as
 * bools indicating whether they are enabled or not.
 *
 * @note The keys map to the lowercase endings of the macros defined in
 * exago_config.h. For example, the key "hiop_sparse" will map to _true_ if
 * the macro EXAGO_ENABLE_HIOP_SPARSE is defined in exago_config.h.
 */
extern const std::unordered_map<std::string, bool> &ExaGOGetBoolConfigOptions();
extern const std::unordered_map<std::string, std::string> &
ExaGOGetStringConfigOptions();

#endif
