/**
 * @file scopflowimpl.h
 * @brief Private header file that defines data and structures for Security constrained Optimal Power Flow application
 */
#ifndef SCOPFLOWIPOPTIMPL_H
#define SCOPFLOWIPOPTIMPL_H

#include <ps.h>
#include <private/psimpl.h>
#include <scopflow.h>
#if defined(PSAPPS_HAVE_IPOPT)
#include <IpStdCInterface.h>
#endif
#include <private/contingencylist.h>
/**
  We want to solve the 2-stage optimization problem with Ns scenarios
  min \sum_{i=0}^Ns f(xi)
  s.t.
     g(xi) = 0              i = 0,Ns  (Equality constraints)
   h^- <= h(xi) <= h^+      i = 0,Ns  (Inequality constraints)
    x^- <= xi <= x+         i = 0,Ns  (Variable bounds)
 -dx <= xi - x0 <= dx       i = 1,Ns  (Coupling constraints between first stage and each scenario)

This can be decomposed for first stage and each scenario as
First stage:
  min f(x0)
  s.t.
     g(x0) = 0              (Equality constraints)
   h^- <= h(x0) <= h^+      (Inequality constraints)
    x^- <= x0 <= x+         (Variable bounds)

For each scenario i \in Ns:
  min f(xi)
  s.t.
     g(xi) = 0              (Equality constraints)
   h^- <= h(xi) <= h^+      (Inequality constraints)
    x^- <= xi <= x+         (Variable bounds)
 -dx <= xi - x0 <= dx       (Coupling constraints between first stage and each scenario)

*/

#define MAX_CONTINGENCIES 100

struct CallBackData{
  void*    prob;
  PetscInt row_node_id;
  PetscInt col_node_id;
  PetscInt typeflag;
};

typedef struct CallBackData *CallBackDataPtr;

struct CCMatrix {
  PetscInt    *rowidx;
  PetscInt    *colptr;
  PetscScalar *values;
};

typedef struct CCMatrix *CCMatrixArray;


 /**
  * @brief private struct for security optimal power flow application
  */
struct _p_SCOPFLOW{
  /* Sizes */
  PetscInt Ns; /* Number of scenarios */
  PetscInt Nx; /* Number of variables */
  PetscInt Nc; /* Number of constraints */
  
  COMM comm; /**< Communicator context */
  OPFLOW *opflows; /* Array of optimal power flow application objects.
  		      Each processor creates ns objects, one for each 
		      scenario */

  PetscBool setupcalled; /* SCOPFLOWSetUp called? */

  PetscBool converged; /* Convergence status */

  char netfile[100]; /* Network data file */

  /* Note that with the above coupling constraints, the Jacobian
     for the coupling constraints is a constant matrix, so we 
     can use just one matrix, instead of one per scenario
  */
  Mat Jcoup;  /* Jacobian for the coupling constraints */
  Mat JcoupT; /* Transpose of the coupling Jacobian */

  PetscBool iscoupling; /* Is each scenario coupled with base scenario? */
  PetscBool first_stage_gen_cost_only; /* Only include the gen cost for first stage only */
  PetscBool ignore_line_flow_constraints; /* Ignore line flow constraints */

#if defined(PSAPPS_HAVE_IPOPT)
  /* IPOPT specific terms */
  IpoptProblem nlp_ipopt; /**< IPOPT problem object */
  enum ApplicationReturnStatus solve_status;
#endif

  PetscInt nnz_jac_g;
  PetscInt nnz_hes;

  Vec lambda_xl;
  Vec lambda_xu;

  PetscScalar obj;

  Vec X,Xl,Xu;
  Vec G,Gl,Gu;
  Vec Lambda;

  PetscInt *xstart; /* Starting location for variables for the scenario in the big X vector */
  PetscInt *gstart; /* Starting location for constraints for the scenario in the big G vector */

  PetscInt *Nxi; /* Number of variables for each scenario */
  PetscInt *Nci; /* Number of constraints for each scenario */

  PetscInt *e_nz_jac_self; /* Number of nonzeros in equality constraint Jacobian for each scenario (Self only) */
  PetscInt *i_nz_jac_self; /* Number of nonzeros in inequality constraint Jacobian for each scenario (self only)*/

  CCMatrixArray e_jac_self; /* Jacobian elements for equality constraint Jacobian for each scenario (self only) */
  CCMatrixArray i_jac_self; /* Jacobian elements for inequality constrained Jacobian for each scenario (self only) */

  PetscInt *e_nz_jac_coupled; /* Number of nonzeros in equality constraint Jacobian for coupling of each scenario with base */
  PetscInt *i_nz_jac_coupled; /* Number of nonzeros in inequality constraint Jacobian for coupling of each scenario with base */

  CCMatrixArray e_jac_coupled; /* Jacobian elements for coupling constraint Jacobian for each scenario (self only) */
  CCMatrixArray i_jac_coupled; /* Jacobian elements for coupling constraint Jacobian for each scenario (self only) */

  PetscInt *nz_hess_self; /* Number of nonzeros in the hessian block for the scenario */

  struct CCMatrix *hess_self; /* Hessian block for each scenario */

  PetscScalar obj_factor; /* The objective factor IPOPT uses in the Hessian evaluation */

  ContingencyList ctgclist;      /* List of contingencies */
  PetscBool       ctgcfileset;   /* Is the contingency file set ? */
  char            ctgcfile[100]; /* Contingency file */

};


#endif
