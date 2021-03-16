#ifndef _SELFCHECK_H
#define _SELFCHECK_H

#include <tcopflow.h>

#define SELFCHECK_NETWORK_CASE9 "case9mod.m"
#define SELFCHECK_NETWORK_CASE118 "case118.m"
#define SELFCHECK_NETWORK_CASE200 "case_ACTIVSg200.m"

#define SELFCHECK_CONTINGENCY_CASE9 "case9.cont"
#define SELFCHECK_SCENARIO_CASE9 "scenarios_9bus.csv"
#define SELFCHECK_SCENARIO_CASE200 "scenarios_200bus.csv"
#define SELFCHECK_NETWORK_CASE9_WIND "case9mod_gen3_wind.m"
#define SELFCHECK_PLOAD "load_P.csv"
#define SELFCHECK_QLOAD "load_Q.csv"

// -- TCOPFLOW Selfcheck --

typedef struct
{
  /// Options passed to the TCOPFLOW driver
  char networkname[PETSC_MAX_PATH_LEN];
  char windgenname[PETSC_MAX_PATH_LEN];
  char qloadprofile[PETSC_MAX_PATH_LEN];
  char ploadprofile[PETSC_MAX_PATH_LEN];
  PetscBool iscoupling;
  PetscReal dt;
  PetscReal duration;
  PetscBool lineflow_constraints;

  /// Attributes of tcopflow solution checked
  PetscInt numiter;
  PetscReal objective;
} ExaGOSelfcheckTCOPFLOWAnswer;

/** Returns 0 if tcopflow conforms to the answer, >0 otherwise. */
extern PetscBool ExaGOSelfcheckTCOPFLOW(TCOPFLOW);
#endif
