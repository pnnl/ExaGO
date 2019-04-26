#ifndef DYNTGOV1_H
#define DYNTGOV1_H

/* Struct declaration for TGOV1 model */
#include <dyn.h>

/* Number of variables for the TGOV1 exciter model */
#define DYNTGOV1_nvar 2

#define TGOV1 "'TGOV1'"

typedef struct _p_DYNTGOV1 *DYNTGOV1;

struct _p_DYNTGOV1{
  PetscInt bus_i; /* Bus number */
  char     id[2]; /* Generator Id */
  PetscScalar T1; /* Valve time constant */
  PetscScalar VMIN; /* Minimum valve position */
  PetscScalar VMAX; /* Maximum valve position */
  PetscScalar T2;   /* Time constant */
  PetscScalar T3;   /* Time constant */
  PetscScalar DT;    /* Speed feedforward gain */
  PetscScalar R;    /* Speed droop constant */
  PetscScalar Pref;  /* Reference power */
  
  PetscInt    X2atmax; /* Flag indicating max. valve position */
  PetscInt    X2atmin; /* Flag indicating min. valve position. */
};

#endif
