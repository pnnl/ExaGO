#ifndef DYNSTAB1_H
#define DYNSTAB1_H

/* Struct declaration for STAB1 model */
#include <dyn.h>

/* Number of variables for the STAB1 exciter model */
#define DYNSTAB1_nvar 3

#define STAB1 "'STAB1'"

typedef struct _p_DYNSTAB1 *DYNSTAB1;

struct _p_DYNSTAB1{
  PetscInt bus_i; /* Bus number */
  char     id[2]; /* Generator Id */
  PetscScalar K_T; /* K/T (sec^-1) */
  PetscScalar T; /* Time constant for washout block (>0)*/
  PetscScalar T1_T3; /* T1/T3 */
  PetscScalar T3; /* Time constant for first lead-lag block (>0) */
  PetscScalar T2_T4; /* T2/T4 */
  PetscScalar T4; /* Time constant for second lead-lag block (>0) */
  PetscScalar Hlim; /* Limit on VOTHSG */
  PetscInt    VOTHSGatmax; /* Flag indicating max. VOTHSG */
  PetscInt    VOTHSGatmin; /* Flag indicating min. VOTHSG */
};

#endif
