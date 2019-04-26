#ifndef DYNCV_H
#define DYNCV_H

#include <dyn.h>

/* Number of variables for the cv model */
#define DYNCv_nvar 2

/* CV model string identifier */
#define CV "'CV'"

typedef struct _p_DYNCv *DYNCv;

struct _p_DYNCv{
  PetscInt bus_i; /* Bus number */
  char id[3];    /* Generator Id  */
  PetscScalar VDs; /* Set point for real-part of voltage */
  PetscScalar VQs; /* Set point for imaginary-part of voltage */
};

#endif
