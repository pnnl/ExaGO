/**
 * @file pflowimpl.h
 * @brief Private header file that defines data and structures for power flow
 * application
 */
#ifndef PFLOWIMPL_H
#define PFLOWIMPL_H

#include <pflow.h>
#include <ps.h>

extern const char *const PFLOWOutputFormatTypes[];

/**
 * @brief Private struct for power flow application
 */
struct _p_PFLOW {
  COMM comm; /**< Communicator context */
  PS ps;     /**< Power system context */

  Vec X; /**< Solution vector */

  Mat Jac; /**< Jacobian */

  SNES snes; /**< Nonlinear solver */

  PetscBool setupcalled; /* PFSetUp called? */
  PetscBool converged;

  OutputFormat outputformat; /* Format for output data */

  PetscBool solutiontops; /* Has PS been updated with solution ? */

  PetscBool
      split_gen_power_within_limits; /* Splits the generator powers (after a
                                      * power flow is solved) by respecing the
                                      * real and reactive power generator
                                      * limits. If 0, it simply splits the
                                      * powers in proportion to the MVAbase
                                      */
};

#endif
