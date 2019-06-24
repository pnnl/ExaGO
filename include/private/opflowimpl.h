/**
 * @file opflowimpl.h
 * @brief Private header file that defines data and structures for Optimal Power Flow application
 */
#ifndef OPFLOWIMPL_H
#define OPFLOWIMPL_H

#include <ps.h>
#include <opflow.h>
#if defined(PSAPPS_HAVE_IPOPT)
#include <IpStdCInterface.h>
#endif
 /** 
  * @brief private struct for optimal power flow application 
  */
struct _p_OPFLOW{
  /* Common fields */
  COMM comm; /**< Communicator context */
  PS   ps;   /**< Power system context */

  Vec  X;    /**< Solution vector */
  Vec  G; /**< Inequality and equality constraint function */

  Vec Xl; /**< Lower bound on solution */
  Vec Xu; /**< Upper bound on solution */

  PetscScalar obj; /**< Objective function */
  Vec gradobj; /**< Gradient of the objective function */

  PetscBool setupcalled; /**< OPFLOWSetUp called? */

  PetscInt Nconeq; /**< Number of equality constraints */
  PetscInt Nconineq; /**< Number of inequality constraints */
  PetscInt Ncon;     /* Total number of constraints (equality + inequality) */
  PetscInt Nvar;     /* Total number of variables */

  PetscInt n; /**< Number of variables */
  PetscInt m; /**< Number of constraints */

  /* For TAO */
  Vec  Ge; /** < Equality constraint function vector */
  Vec  Gi; /** < Inequality constraint function vector */
  
  Mat Jac_Ge; /* Equality constraint Jacobian */
  Mat Jac_Gi; /* Inequality constraint Jacobian */

  Tao nlp;    /* Optimization problem */

  PetscBool converged; // Convergence status

  /* For IPOPT */
  Vec Gl; /**< Lower bound on G */
  Vec Gu; /**< Upper bound on G */

  Mat  Jac;  /* Jacobian of constraints */
  Mat  Hes;  /* Lagrangian Hessian */  
  
  PetscScalar *x; /**< Solution array - same as the array for X */
  
  PetscInt nnz_jac_g; /**< Number of nonzeros in the jacobian of the constraints */
  PetscInt nnz_hes; /**< Number of nonzeros in the Lagrangian Hessian */

  /* Lagrange multipliers */
  Vec lambda_g;
  Vec lambda_xl;
  Vec lambda_xu;

#if defined(PSAPPS_HAVE_IPOPT)
  /* Ipopt specific terms */
  IpoptProblem nlp_ipopt; /**< Ipopt solver */
  enum ApplicationReturnStatus solve_status;
#endif


};
/**
 * @brief Sets the bounds on variables and constraints
 * @param [in] OPFLOW opflow - the OPFLOW object
 * @param [out] Vec Xl     - vector of lower bound on variables
 * @param [out] Vec Xu     - vector of upper bound on variables
 * @param [out] Vec Gl     - vector of lower bound on constraints
 * @param [out] Vec Gu     - vector of lower bound on constraints
 */
extern PetscErrorCode OPFLOWSetVariableandConstraintBounds(OPFLOW,Vec,Vec,Vec,Vec);
/**
 * @brief Sets the initial guess for the optimization
 * @param [in] OPFLOW opflow - the OPFLOW object
 * @param [out] Vec X     - initial guess
 */
extern PetscErrorCode OPFLOWSetInitialGuess(OPFLOW,Vec);
#endif
