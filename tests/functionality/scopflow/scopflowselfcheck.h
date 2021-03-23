#ifndef _SCOPFLOWSELFCHECK_H
#define _SCOPFLOWSELFCHECK_H

#include <scopflow.h>

#define SELFCHECK_NETWORK_CASE9 "case9mod.m"
#define SELFCHECK_NETWORK_CASE118 "case118.m"
#define SELFCHECK_NETWORK_CASE200 "case_ACTIVSg200.m"

// -- SCOPFLOW Selfcheck --

#define SELFCHECK_CONTINGENCY_CASE9 "case9.cont"
#define SELFCHECK_CONTINGENCY_CASE200 "case_ACTIVSg200.cont"
#define SELFCHECK_SCENARIO_CASE9 "scenarios_9bus.csv"
#define SELFCHECK_SCENARIO_CASE200 "scenarios_200bus.csv"
#define SELFCHECK_NETWORK_CASE9_WIND "case9mod_gen3_wind.m"
#define SELFCHECK_PLOAD "load_P.csv"
#define SELFCHECK_QLOAD "load_Q.csv"

typedef struct
{
  /// Options passed to scopflow driver that are:
  ///  - applicable to all scopflow runs
  char solver[PETSC_MAX_PATH_LEN];
  char networkname[PETSC_MAX_PATH_LEN];
  char scenarioname[PETSC_MAX_PATH_LEN];
  char modelinit[PETSC_MAX_PATH_LEN];
  char genbusvoltage[PETSC_MAX_PATH_LEN];
  PetscInt numcontingencies;
  char contingencyname[PETSC_MAX_PATH_LEN];
  PetscInt mode; /** preventive(0) -default ; and corrective(1) */

  PetscBool multiperiod;
  ///  - applicable to multi-period
  char ploadprofile[PETSC_MAX_PATH_LEN];
  char qloadprofile[PETSC_MAX_PATH_LEN];
  PetscReal dt;
  PetscReal duration;

  /// Attributes of scopflow solution checked
  PetscInt numiter;
  PetscReal objective;
} ExaGOSelfcheckSCOPFLOWAnswer;

/** Returns 0 if scopflow conforms to the answer, >0 otherwise. */
extern PetscBool ExaGOSelfcheckSCOPFLOW(SCOPFLOW);
#endif
