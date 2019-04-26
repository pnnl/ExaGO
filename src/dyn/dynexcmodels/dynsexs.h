#ifndef DYNSEXS_H
#define DYNSEXS_H

/* Struct declaration for SEXS model */
#include <dyn.h>

/* Number of variables for the SEXS exciter model */
#define DYNSEXS_nvar 2

#define SEXS "'SEXS'"

typedef struct _p_DYNSEXS *DYNSEXS;

struct _p_DYNSEXS{
  PetscInt bus_i; /* Bus number */
  char     id[3]; /* Generator Id */
  PetscScalar TA_TB; /* Time-constant ratio TA/TB */
  PetscScalar TB; /* Time-constant TB */
  PetscScalar K;  /* Exciter gain */
  PetscScalar TE; /* Exciter time constant */
  PetscScalar Emax; /* Max. allowed E (pu on Efd base) */
  PetscScalar Emin; /* Min. allowed E (pu on Efd base) */
  PetscScalar Vref;  /* Reference voltage */
  
  PetscInt    Eatmax; /* Flag indicating E is at max. */
  PetscInt    Eatmin; /* Flag indicating E is at min. */
};

#endif
