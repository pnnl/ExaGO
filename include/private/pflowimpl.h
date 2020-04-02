/**
 * @file pflowimpl.h
 * @brief Private header file that defines data and structures for power flow application
 */
#ifndef PFLOWIMPL_H
#define PFLOWIMPL_H

#include <ps.h>
#include <pflow.h>

/**
 * @brief Private struct for power flow application
 */
struct _p_PFLOW{
  COMM comm; /**< Communicator context */
  PS   ps;   /**< Power system context */

  Vec  X;    /**< Solution vector */

  Mat  Jac;  /**< Jacobian */

  SNES snes; /**< Nonlinear solver */

  PetscBool setupcalled; /* PFSetUp called? */
  PetscBool converged;

  PetscBool  split_gen_power_within_limits; /* Splits the generator powers (after a power flow is solved) by respecing the real and reactive power generator limits. If 0, it simply splits the powers in proportion to the MVAbase
					     */
};

#endif
