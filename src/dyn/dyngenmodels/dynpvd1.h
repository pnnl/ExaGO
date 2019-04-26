#ifndef DYNPVD1_H
#define DYNPVD1_H

#include <dyn.h>

/* Distributed PV system model
   http://prod.sandia.gov/techlib/access-control.cgi/2013/138876.pdf
   Default values from NERC report draft on DER modeling
   http://www.nerc.com/comm/PC/LoadModelingTaskForceDL/Reliability_Guideline_-_DER_Modeling_Parameters_-_2017-04-04_-_FINAL_DRAFT.pdf
*/

/* PVD1 model string identifier */
#define PVD1 "'PVD1'"

typedef struct _p_DYNPvd1 *DYNPvd1;

struct _p_DYNPvd1{
  PetscInt bus_i; /* Bus number */
  char id[3];    /* Generator Id  */
  PetscInt pqflag; /* Priority to reactive or active current (reactive current = 0, active current = 1) */
  PetscScalar xc; /* Line drop compensation reactance (pu) (0) */
  PetscScalar qmx; /* Maximum reactive power command (pu) */
  PetscScalar qmn; /* Minimum reactive power command (pu) */
  PetscScalar v0; /* Lower limit of deadband for voltage droop response */
  PetscScalar v1; /* Uppoer limit of deadband for voltage droop response */
  PetscScalar dqdv; /* Voltage droop characteristic */
  PetscScalar fdbd; /* Overfrequency deadband for governor response (pu) */
  PetscScalar ddn;  /* Down regulation droop gain (pu) */
  PetscScalar imax; /* Apparent current limit (pu) */
  PetscScalar vt0; /* Voltage tripping response curve point 0 (pu) */
  PetscScalar vt1; /* Voltage tripping response curve point 1 (pu) */
  PetscScalar vt2; /* Voltage tripping response curve point 2 (pu) */
  PetscScalar vt3; /* Voltage tripping response curve point 3 (pu) */
  PetscInt    vrflag; /* Voltage tripping method */
  PetscScalar ft0; /* Frequency tripping response curve point 0 (Hz) */
  PetscScalar ft1; /* Frequency tripping response curve point 1 (Hz) */
  PetscScalar ft2; /* Frequency tripping response curve point 2 (Hz) */
  PetscScalar ft3; /* Frequency tripping response curve point 3 (Hz) */
  PetscInt    frflag; /* Frequency tripping method */
  PetscScalar tg; /* Inverter current lag time constant (sec) */
  PetscScalar tf; /* Frequency transducer time constant (sec) */
  PetscScalar vtmax; /* Voltage limit used in high voltage reactive power logic (pu) */
  PetscScalar Ivpnt1; /* High voltage point for low voltage active current management function (pu) */
  PetscScalar Ivpnt0; /* Low voltage point for low voltage active current management function (pu) */
  PetscScalar qmin;   /* Limit in high voltage reactive power (pu) */
  PetscScalar Khv; /* Acceleration factor used in high voltage reactive power logic (pu) */

  PetscScalar Pref; /* Real power reference (P from power flow) */
  PetscScalar Qref; /* Reactive power reference (Q from power flow) */
  PetscScalar fvl; /* low voltage tripping fraction */
  PetscScalar fvh; /* High voltage tripping fraction */
};

#endif
