#ifndef _OPFLOWSELFCHECK_H
#define _OPFLOWSELFCHECK_H

#include <opflow.h>

#define SELFCHECK_NETWORK_CASE9 "case9mod.m"
#define SELFCHECK_NETWORK_CASE118 "case118.m"
#define SELFCHECK_NETWORK_CASE200 "case_ACTIVSg200.m"

// -- OPFLOW Selfcheck --

/** Description of OPFLOW problem used for -selfcheck flag */
typedef struct
{
  /// Options passed to opflow driver
  char solver[PETSC_MAX_PATH_LEN];
  char model[PETSC_MAX_PATH_LEN];
  char networkname[PETSC_MAX_PATH_LEN];
  bool  isloadloss_active; // Are load loss variables active?
  bool  ispowerimb_active; // Are power imbalance variables active?

  /// Attributes of opflow solution checked
  PetscInt numiter;
  PetscReal objective;
} ExaGOSelfcheckOPFLOWAnswer;

/** Returns 0 if opflow conforms to the answer, >0 otherwise. */
extern PetscErrorCode ExaGOSelfcheckOPFLOW(OPFLOW);

#endif
