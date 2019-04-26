/**
 * @file dynopflowimpl.h
 * @brief Private header file that defines objects for dynamics optimization application
 */
#ifndef DYNOPFLOWIMPL_H
#define DYNOPFLOWIMPL_H

#include <ps.h>
#include <opflow.h>
#include <dyn.h>
#include <dynopflow.h>
#include "IpStdCInterface.h" 

 /** 
  * @brief Private struct for dynamics constrained optimization application
  */
struct _p_DYNOPFLOW{
  COMM comm; /**<  Communicator context */
  OPFLOW opflow; /**<  Optimal power flow context */
  DYN    dyn; /**<  Dynamics application context */
  PetscBool setupcalled; /**<  Is setup called? */
  PetscInt Nconopflow; /**<  Number of OPFLOW constraints */
  PetscInt Ncondyn; /**<  Number of dynamic constraints */
  PetscInt Ncon; /**<  Total number of constraints */

  /* For sensitivity calculation */
  Vec *lambda; /**<  Gradients of cost functions w.r.t initial conditions */
  Vec *mu; /**<  Gradients of cost functions w.r.t parameters */
  Vec *dy0dp; /**<  Partial derivatives of initial conditions w.r.t parameters */
  Mat  Jacp; /**<  Jacobian of the DAE RHS w.r.t parameters */
  Vec  costintegral; /**<  Cost Integral */

  /* Cost function 
     r = scal*f^exp -> where f is the cost function and r is the constraint
     that gets incorporated in the optimization
  */
  PetscReal scal; /**<  Scaling of the cost function */
  PetscReal exp; /**<  Exponent for smoothing cost function */
  PetscReal eta; /**<  Constraint bound (r <= \eta) */

  /* Constraint Jacobian nonzero information */
  PetscInt nnz_opflow_jac_g; /**<  Number of nonzeros in the constraint Jacobian of the OPFLOW part */
  PetscInt nnz_dyn_jac_g; /**<  Number of nonzeros in the constraint Jacobian due to dynamic constraints */
  PetscInt nnz_jac_g; /**<  Total number of nonzeros in the dynopflow constraint jacobian (=nnz_opflow_jac_g + nnz_dyn_jac_g) */

  /* Limits for frequency */
  PetscReal freq_max; /**< Max. frequency */
  PetscReal freq_min; /**< Min. frequency */
  
  PetscBool initsolve; /**<  Initialization flag (solve OPFLOW only if TRUE) */

  PetscBool solvesimultaneous; /**<  Solve simultaneously (1), cutting plane (0) */

  PetscBool sum_freq_cons; /**<  When this flag is set, a single constraint is used for the frequency constraints */

  /* For cutting plane method */
  PetscInt cpmaxit; /**<  Max. iterations for cutting plane method */
  Vec      Xpre; /**<  Previous iterate */
  Vec      Xtemp; /**<  Temporary vector used in operations */
  Vec      Gdynpre; /**<  Dynamic constraints from previous iterate */
};

#endif


