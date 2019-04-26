#ifndef DYNEXST1_H
#define DYNEXST1_H

/* Struct declaration for EXST1 model */
#include <dyn.h>

/* Number of variables for the EXST1 exciter model */
#define DYNEXST1_nvar 4

#define EXST1 "'EXST1'"

typedef struct _p_DYNEXST1 *DYNEXST1;

struct _p_DYNEXST1{
  PetscInt bus_i; /* Bus number */
  char     id[3]; /* Generator Id */

  /* Add fields to hold data for the model */


  PetscScalar TR; /* Transducer time constant */
  PetscScalar KA; /* Amplifier gain */
  PetscScalar Vimin; /* Max. allowed VI */
  PetscScalar Vimax; /* Min. allowed VI */
  PetscScalar TC; /* The lead time constant of the lead-lag compensation block */
  PetscScalar TB; /* The lag time constant of the lead-lag compensation block */
  PetscScalar TA; /* Amplifier time constant */
  PetscScalar TF;    /* Feedback stabilizer time constant */
  PetscScalar KF;    /* Feedback stabilizer gain */
  PetscScalar VRmax; /* Max. allowed VR */
  PetscScalar VRmin; /* Min. allowed VR */
  PetscScalar KC;    /* Commuting reactance */
  PetscScalar Vref;  /* Reference voltage */
  
  PetscInt    VIatmax; /* Flag indicating VI is at max. */
  PetscInt    VIatmin; /* Flag indicating VI is at min. */
  PetscInt    Efdatmax; /* Flag indicating Efd is at max. */
  PetscInt    Efdatmin; /* Flag indicating Efd is at min. */
  PetscScalar Efdmin; /* Minimum allowed Efd THIS MIGHT NOT BE NECESSARY AS A PARAMETER*/
  PetscScalar Efdmax; /* Maximum allowed Efd THIS MIGHT NOT BE NECESSARY AS A PARAMETER*/

};

extern PetscErrorCode DYNExcModelCreate_EXST1(DYNExcModel);
#endif
