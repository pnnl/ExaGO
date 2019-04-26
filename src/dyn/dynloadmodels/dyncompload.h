#ifndef DYNCOMPLOAD_H
#define DYNCOMPLOAD_H

/* TO DO
  -- 1: Include the impact of MBase and Pmult
*/

/* Comploadosite load model comprising of ZIP load and motor */
#include <dyn.h>

/* COMPLOAD model string identifier */
#define COMPLOAD "'COMPLOAD'"

typedef struct _p_DYNCOMPLOAD *DYNCOMPLOAD;

/*
  BUSI "COMPLOAD" ID PL QL IP IQ YP YQ IT RA XA XM R1 X1 R2 X2 E1 SE1 E2 SE2 MBASE PMULT \
  H VI TI TB D Tnom
*/
struct _p_DYNCOMPLOAD{
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
  PetscScalar pmot; /* Fraction of motor real power consumption */
  PetscBool   motstatus; /* Status of motor model, OFF if pl+ip+yp = 1 */
  PetscBool   cages; /* Single or double cage */
  /* Motor model parameters */
  PetscInt IT;
  PetscScalar Ra;
  PetscScalar Xa;
  PetscScalar Xm;
  PetscScalar R1;
  PetscScalar X1;
  PetscScalar R2;
  PetscScalar X2;
  PetscScalar E1;
  PetscScalar SE1;
  PetscScalar E2;
  PetscScalar SE2;
  PetscScalar Mbase;
  PetscScalar Pmult;
  PetscScalar H;
  PetscScalar Vi;
  PetscScalar Ti;
  PetscScalar Tb;
  PetscScalar D;
  PetscScalar Tnom;
  /* Motor model auxillary parameters */
  PetscScalar s0; /* Initial slip -- This is currently set to 0.05 during initiliaziation
		     for convenience only. Later, the code should be updated to 
		     calculate slip.
		  */
  PetscScalar Bcomp; /*  Shunt susceptance to compensate the difference between
		         reactive power load specified by power flow and 
		         reactive power drawn by the motor */
  PetscScalar X0;
  PetscScalar Xp;
  PetscScalar Tp0;
  PetscScalar Tmech;

  /* Aux. variables used in motor tripping event detection */
  PetscScalar t_timer_start;    /* Time at which the time starts */
  PetscInt    low_voltage_flag; /* Flag to indicate low voltage V < VI */
  PetscInt    mot_on;   /* Flag to indicate motor trip (0 = OFF, 1 = ON) */
};

#endif
