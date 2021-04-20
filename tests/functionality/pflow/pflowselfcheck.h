#ifndef _PFLOWSELFCHECK_H
#define _PFLOWSELFCHECK_H

#include <pflow.h>

#define SELFCHECK_NETWORK_CASE9 "case9mod.m"
#define SELFCHECK_NETWORK_CASE118 "case118.m"
#define SELFCHECK_NETWORK_CASE200 "case_ACTIVSg200.m"

// -- PFLOW Selfcheck --

/** Description of PFLOW problem used for -selfcheck flag */
typedef struct
{
  /// Options passed to the pflow driver
  char networkname[PETSC_MAX_PATH_LEN];

  /// Attributes of pflow solution checked
  PetscInt numiter;
} ExaGOSelfcheckPFLOWAnswer;

/** Returns 0 if opflow conforms to the answer, >0 otherwise. */
extern PetscErrorCode ExaGOSelfcheckPFLOW(PFLOW);

#endif
