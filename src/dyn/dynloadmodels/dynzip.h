#ifndef DYNZIP_H
#define DYNZIP_H

#include <dyn.h>

/* Number of variables for the zip model */
#define DYNZIP_nvar 0

/* ZIP model string identifier */
#define ZIP "'ZIP'"

typedef struct _p_DYNZIP *DYNZIP;

struct _p_DYNZIP{
  PetscInt bus_i;  /* Bus number */
  char id[3];      /* Load Id  */
  PetscScalar Vm0; /* Voltage magnitude at simulation start */
  PetscScalar pl;  /* Real part of constant power portion */
  PetscScalar ql;  /* Imaginary part of constant power portion */
  PetscScalar ip;  /* Real part of constant current portion */
  PetscScalar iq;  /* Real part of constant current portion */
  PetscScalar yp;  /* Real part of constant impedance portion */
  PetscScalar yq;  /* Real part of constant impedance portion */
  PetscScalar Vm_thresh; /* Voltage magnitude threshold below which all constant power and constant current loads are converted to constant impedance */
};

#endif
