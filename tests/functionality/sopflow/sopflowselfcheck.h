#ifndef _SOPFLOWSELFCHECK_H
#define _SOPFLOWSELFCHECK_H

#include <sopflow.h>

#define SELFCHECK_NETWORK_CASE9 "case9mod.m"
#define SELFCHECK_NETWORK_CASE118 "case118.m"
#define SELFCHECK_NETWORK_CASE200 "case_ACTIVSg200.m"

// -- SOPFLOW Selfcheck --

#define SELFCHECK_CONTINGENCY_CASE9 "case9.cont"
#define SELFCHECK_CONTINGENCY_CASE200 "case_ACTIVSg200.cont"
#define SELFCHECK_SCENARIO_CASE9 "scenarios_9bus.csv"
#define SELFCHECK_SCENARIO_CASE200 "scenarios_200bus.csv"
#define SELFCHECK_NETWORK_CASE9_WIND "case9mod_gen3_wind.m"
#define SELFCHECK_PLOAD "load_P.csv"
#define SELFCHECK_QLOAD "load_Q.csv"

typedef struct
{
  /// Options passed to sopflow driver that are:
  ///  - applicable to all sopflow runs
  char solver[PETSC_MAX_PATH_LEN];
  char networkname[PETSC_MAX_PATH_LEN];
  char scenarioname[PETSC_MAX_PATH_LEN];
  char modelinit[PETSC_MAX_PATH_LEN];
  char genbusvoltage[PETSC_MAX_PATH_LEN];
  PetscInt numscenarios;

  PetscBool multicontingency;
  ///  - applicable to muti-contingency
  PetscInt numcontingencies;
  char contingencyname[PETSC_MAX_PATH_LEN];

  PetscBool multiperiod;
  ///  - applicable to multi-period
  char ploadprofile[PETSC_MAX_PATH_LEN];
  char qloadprofile[PETSC_MAX_PATH_LEN];
  PetscReal dt;
  PetscReal duration;

  /// Attributes of sopflow solution checked
  PetscInt numiter;
  PetscReal objective;
} ExaGOSelfcheckSOPFLOWAnswer;

/** Returns 0 if sopflow conforms to the answer, >0 otherwise. */
extern PetscErrorCode ExaGOSelfcheckSOPFLOW(SOPFLOW);
#endif
