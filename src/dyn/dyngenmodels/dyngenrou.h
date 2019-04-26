#ifndef DYNGENROU_H
#define DYNGENROU_H

#include <dyn.h>

/* Number of variables for the genrou model */
#define DYNGenrou_nvar 8

/* GENROU model string identifier */
#define GENROU "'GENROU'"

typedef struct _p_DYNGenrou *DYNGenrou;

struct _p_DYNGenrou{
  PetscInt bus_i; /* Bus number */
  char id[3];    /* Generator Id  */
  PetscScalar Td0p;  /* d-axis open circuit transient time constant */
  PetscScalar Td0dp; /* d-axis open circuit sub-transient time constant */
  PetscScalar Tq0p;  /* q-axis open circuit transient time constant */
  PetscScalar Tq0dp; /* q-axis open circuit sub-transient time constant */ 
  PetscScalar H;      /* Inertial constant */
  PetscScalar D;      /* Damping coefficient */
  PetscScalar Xd;     /* d-axis reactance */
  PetscScalar Xdp;    /* d-axis transient reactance */
  PetscScalar Xq;     /* q-axis reactance */
  PetscScalar Xqp;    /* q-axis transient reactance */
  PetscScalar Xddp;   /* d-axis sub-transient reactance = Xqdp */
  PetscScalar Xl;     /* Leakage reactance */
  PetscScalar S1;     /* Saturation constant */
  PetscScalar S2;     /* Saturation constant */
  PetscScalar Pm;     /* mechanical power input */
  PetscScalar Efd;    /* field voltage */
  PetscInt    freqatmax; /* Flag to indicate freq at max */
  PetscInt    freqatmin; /* Flag to indicate freq at min */
};

#endif
