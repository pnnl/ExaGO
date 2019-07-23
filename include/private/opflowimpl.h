/**
 * @file opflowimpl.h
 * @brief Private header file that defines data and structures for Optimal Power Flow application
 */
#ifndef OPFLOWIMPL_H
#define OPFLOWIMPL_H

#include <ps.h>
#include <private/psimpl.h>
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

  Vec  X,localX;    /* Global and local solution vector */
  Vec  G; /**< Inequality and equality constraint function */

  Vec Xl; /**< Lower bound on solution */
  Vec Xu; /**< Upper bound on solution */

  PetscScalar obj; /**< Objective function */
  Vec gradobj; /**< Gradient of the objective function */

  PetscBool setupcalled; /* OPFLOWSetUp called? */

  PetscInt nconeq, Nconeq;     /* Local and global number of equality constraints, excluding ghosts! */
  PetscInt nconineq, Nconineq; /* Local and global number of inequality constraints */
  PetscInt Ncon;               /* Total number of constraints (equality + inequality) */
  PetscInt nvar,Nvar;          /* Total number of local and global variables, excluding ghosts! */

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
/*
  * @breif Evalulates the equality constraints for the optimal power flow
  * @param [in] Tao - the Tao nlp solver object
  * @param [in] X   - the current iterate
  * @param [in] ctx - application data set with TaoSetEqualityConstraintsRoutine
  * @param [out] Ge  - vector of equality constraints
*/
extern PetscErrorCode OPFLOWEqualityConstraintsFunction(Tao nlp,Vec X,Vec Ge,void* ctx);
/*
  * @breif Returns a distributed Jacobian for the equality constraints
  * @param [in] opflow - the optimal power flow application object
  * @param [out] mat - the jacobian of the equality constraints
*/
extern PetscErrorCode OPFLOWCreateEqualityConstraintsJacobian(OPFLOW opflow,Mat *mat);
/*
  * @breif  Sets the nonzero values for the equality constraints Jacobian
  * @param [in] nlp - Tao nlp solver object
  * @param [in] X   - the current iterate
  * @param [in] ctx - application data set with OPFLOWEqualityConstraintsJacobianRoutine
  * @param [out] Je - Jacobian of equality constraints
  * @param [out] Je_pre - Preconditioner for equality constraints
*/
extern PetscErrorCode OPFLOWEqualityConstraintsJacobianFunction(Tao nlp,Vec X,Mat Je, Mat Je_pre, void* ctx);
/*
  * @breif Evalulates the inequality constraints for the optimal power flow
  * @param [in] Tao - the Tao nlp solver object
  * @param [in] X   - the current iterate
  * @param [in] ctx - application data set with TaoSetEqualityConstraintsRoutine
  * @param [out] Gi  - vector of equality constraints
*/
extern PetscErrorCode OPFLOWInequalityConstraintsFunction(Tao nlp,Vec X,Vec Gi,void* ctx);
/*
  * @breif Returns a distributed Jacobian for the inequality constraints
  * @param [in] opflow - the optimal power flow application object
  * @param [out] mat - the jacobian of the inequality constraints
*/
extern PetscErrorCode OPFLOWCreateInequalityConstraintsJacobian(OPFLOW opflow,Mat *mat);
/*
  * @breif  Sets the nonzero values for the equality constraints Jacobian
  * @param [in] nlp - Tao nlp solver object
  * @param [in] X   - the current iterate
  * @param [in] ctx - application data set with OPFLOWEqualityConstraintsJacobianRoutine
  * @param [out] Ji - Jacobian of inequality constraints
  * @param [out] Ji_pre - Preconditioner for inequality constraints
*/
extern PetscErrorCode OPFLOWInequalityConstraintsJacobianFunction(Tao nlp, Vec X, Mat Ji, Mat Ji_pre, void* ctx);

#endif
