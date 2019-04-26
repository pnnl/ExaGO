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
};

#endif
