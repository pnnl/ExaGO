/**
 * @file scopflowimpl.h
 * @brief Private header file that defines data and structures for Security constrained Optimal Power Flow application
 */
#ifndef SCOPFLOWIMPL_H
#define SCOPFLOWIMPL_H

#include <ps.h>
#include <private/psimpl.h>
#include <opflow.h>
#include <scopflow.h>
#include <private/contingencylist.h>

#define SCOPFLOWSOLVERSMAX 3

struct _p_SCOPFLOWSolverList{
  char name[32]; /* Name of the solver */
  PetscErrorCode (*create)(SCOPFLOW); /* Solver object creation routine */
};

struct _p_SCOPFLOWSolverOps {
  PetscErrorCode (*destroy)(SCOPFLOW);
  PetscErrorCode (*setup)(SCOPFLOW);
  PetscErrorCode (*solve)(SCOPFLOW);
  PetscErrorCode (*getobjective)(SCOPFLOW,PetscReal*);
  PetscErrorCode (*getbasecasesolution)(SCOPFLOW,Vec*);
  PetscErrorCode (*getconvergencestatus)(SCOPFLOW,PetscBool*);
  PetscErrorCode (*getconstraints)(SCOPFLOW,Vec*);
  PetscErrorCode (*getconstraintmultipliers)(SCOPFLOW,Vec*);
};


/**
 * @brief private struct for security optimal power flow application
*/
struct _p_SCOPFLOW{
  /* Sizes */
  PetscInt ns,Ns; /* Number of local and global (total) scenarios */
  PetscInt nx,Nx; /* Local and global (total) number of variables */
  PetscInt Ncon,ncon; /* Number of constraints */
  PetscInt Nconeq,nconeq; /* Local and global number of equality constraints */
  PetscInt nconineq, Nconineq; /* Local and global number of inequality constraints */
  PetscInt *nconineqcoup;      /* Number of inequality coupling constraints */
  
  COMM comm; /**< Communicator context */
  OPFLOW *opflows; /* Array of optimal power flow application objects.
  		      Each processor creates ns objects, one for each 
		      scenario */

  PetscBool setupcalled; /* SCOPFLOWSetUp called? */

  PetscBool converged; /* Convergence status */

  char netfile[100]; /* Network data file */

  Vec  X,localX;    /* Global and local solution vector */
  Vec  G; /**< Inequality and equality constraint function */
  Vec  Ge,Gelocal; /** < Equality constraint function vector (global and local) */
  Vec  Gi; /** < Inequality constraint function vector */
  Vec  Lambdai; 
  Vec  Lambdae,Lambdaelocal;
  Vec  Lambda;

  Mat  Jac; /* Jacobian for equality and inequality constraints */
  Mat  Jac_Ge; /* Equality constraint Jacobian */
  Mat  Jac_Gi; /* Inequality constraint Jacobian */

  Mat  Hes;  /* Lagrangian Hessian */
  
  Vec Xl; /**< Lower bound on solution */
  Vec Xu; /**< Upper bound on solution */

  Vec Gl; /**< Lower bound on G */
  Vec Gu; /**< Upper bound on G */

  PetscScalar obj; /**< Objective function */
  Vec gradobj; /**< Gradient of the objective function */
  
  PetscScalar obj_factor; /* IPOPT scales the objective hessian part with this factor. For all other solvers, unless it is set, obj_factor = 1.0. */

  Mat *Jcoup;  /* Jacobian for the coupling constraints (one per scenario) */
  Mat *JcoupT; /* Transpose of the coupling Jacobian (one per scenario)*/

  PetscBool iscoupling; /* Is each scenario coupled with base scenario? */
  PetscBool first_stage_gen_cost_only; /* Only include the gen cost for first stage only */
  PetscBool replicate_basecase; /* Replicate base case for all scenarios */

  void *solver; /* Solver object */
  struct _p_SCOPFLOWSolverOps solverops;

  /* List of solvers registered */
  struct _p_SCOPFLOWSolverList SCOPFLOWSolverList[SCOPFLOWSOLVERSMAX];
  PetscInt nsolversregistered;
  PetscBool SCOPFLOWSolverRegisterAllCalled;

  void    *solverdata;

  /* Data for contingencies */
  ContingencyList ctgclist;      /* List of contingencies */
  PetscBool       ctgcfileset;   /* Is the contingency file set ? */
  char            ctgcfile[100]; /* Contingency file */

};

/* Register all SCOPFLOW solvers */
extern PetscErrorCode SCOPFLOWSolverRegisterAll(SCOPFLOW);
extern PetscErrorCode SCOPFLOWGetConstraints(SCOPFLOW,Vec*);
extern PetscErrorCode SCOPFLOWGetConstraintMultipliers(SCOPFLOW,Vec*);


#endif
