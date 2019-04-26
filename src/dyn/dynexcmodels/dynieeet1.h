#ifndef DYNIEEET1_H
#define DYNIEEET1_H

/* Struct declaration for IEEET1 model */
#include <dyn.h>

/* Number of variables for the IEEET1 exciter model */
#define DYNIEEET1_nvar 4

#define IEEET1 "'IEEET1'"

typedef struct _p_DYNIEEET1 *DYNIEEET1;

struct _p_DYNIEEET1{
  PetscInt bus_i; /* Bus number */
  char     id[3]; /* Generator Id */
  PetscScalar TR; /* Transducer time constant */
  PetscScalar KA; /* Amplifier gain */
  PetscScalar TA; /* Amplifier time constant */
  PetscScalar VRmax; /* Max. allowed VR */
  PetscScalar VRmin; /* Min. allowed VR */
  PetscScalar KE;    /* Exciter gain */
  PetscScalar TE;    /* Exciter time constant */
  PetscScalar KF;    /* Feedback stabilizer gain */
  PetscScalar TF;    /* Feedback stabilizer time constant */
  PetscScalar E1;    /* Efd1 on Efd-SE(Efd) curve */
  PetscScalar SE1;   /* SE1 on Efd-SE(Efd) curve */
  PetscScalar E2;    /* Efd2 on Efd-SE(Efd) curve */
  PetscScalar SE2;   /* SE2 on Efd-SE(Efd) curve */
  PetscScalar Vref;  /* Reference voltage */
  
  PetscScalar satA; /* Saturation constant */
  PetscScalar satB; /* Saturation constant */
  PetscScalar Efdthresh; /* Threshold for Efd saturation */

  PetscInt    VRatmax; /* Flag indicating VR is at max. */
  PetscInt    VRatmin; /* Flag indicating VR is at min. */
};

#endif
