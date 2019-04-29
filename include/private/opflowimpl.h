/**
 * @file opflowimpl.h
 * @brief Private header file that defines data and structures for Optimal Power Flow application
 */
#ifndef OPFLOWIMPL_H
#define OPFLOWIMPL_H

#include <ps.h>
#include <opflow.h>

 /** 
  * @brief private struct for optimal power flow application 
  */
struct _p_OPFLOW{
  COMM comm; /**< Communicator context */
  PS   ps;   /**< Power system context */

  Vec  X;    /**< Solution vector */
  PetscScalar *x; /**< Solution array - same as the array for X */

  Vec  G; /**< Inequality and equality constraint function */
  Vec Gl; /**< Lower bound on G */
  Vec Gu; /**< Upper bound on G */

  Vec  Ge; /** < Equality constraint function vector */
  Vec  Gi; /** < Inequality constraint function vector */

  Vec Gel; /**< Lower bound on Ge */
  Vec Geu; /**< Upper bound on Ge */
  Vec Gil; /**< Lower bound on Gi */
  Vec Giu; /**< Upper bound on Gi */

  PetscScalar obj; /**< Objective function */

  Vec gradobj; /**< Gradient of the objective function */

  Vec Xl; /**< Lower bound on solution */
  Vec Xu; /**< Upper bound on solution */

  Mat  Jac;  /* Jacobian of constraints */

  Mat Jac_Ge; /* Equality constraint Jacobian */
  Mat Jac_Gi; /* Inequality constraint Jacobian */

  Mat  Hes;  /* Lagrangian Hessian */  

  PetscInt Nconeq; /**< Number of equality constraints */
  PetscInt Nconineq; /**< Number of inequality constraints */
  PetscInt Ncon;     /* Total number of constraints (equality + inequality) */
  PetscInt Nvar;     /* Total number of variables */

  PetscInt n; /**< Number of variables */
  PetscInt m; /**< Number of constraints */
  PetscInt nnz_jac_g; /**< Number of nonzeros in the jacobian of the constraints */
  PetscInt nnz_hes; /**< Number of nonzeros in the Lagrangian Hessian */

  PetscBool setupcalled; /**< OPFLOWSetUp called? */

  /* Lagrange multipliers */
  Vec lambda_g;
  Vec lambda_xl;
  Vec lambda_xu;

  /* Ipopt specific terms */
  //  IpoptProblem nlp; /**< Ipopt solver */
  Tao nlp;

  PetscBool converged; // Convergence status
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
