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

#define SCOPFLOWSOLVERSMAX 10

struct _p_SCOPFLOWSolverList{
  char name[32]; /* Name of the solver */
  PetscErrorCode (*create)(SCOPFLOW); /* Solver object creation routine */
};

struct _p_SCOPFLOWSolverOps {
  PetscErrorCode (*destroy)(SCOPFLOW);
  PetscErrorCode (*setup)(SCOPFLOW);
  PetscErrorCode (*solve)(SCOPFLOW);
  PetscErrorCode (*getobjective)(SCOPFLOW,PetscReal*);
  PetscErrorCode (*getsolution)(SCOPFLOW,PetscInt,Vec*);
  PetscErrorCode (*getconvergencestatus)(SCOPFLOW,PetscBool*);
  PetscErrorCode (*getconstraints)(SCOPFLOW,PetscInt,Vec*);
  PetscErrorCode (*getconstraintmultipliers)(SCOPFLOW,PetscInt,Vec*);
};


/**
 * @brief private struct for security optimal power flow application
*/
struct _p_SCOPFLOW{
  /* Sizes */
  PetscInt nc,Nc; /* Number of local and global (total) contingencies */
  PetscInt nx,Nx; /* Local and global (total) number of variables */
  PetscInt Ncon,ncon; /* Number of constraints */
  PetscInt Nconeq,nconeq; /* Local and global number of equality constraints */
  PetscInt nconineq, Nconineq; /* Local and global number of inequality constraints */
  PetscInt *nconineqcoup;      /* Number of inequality coupling constraints */
  PetscInt Nconcoup;           /* Number of coupling constraints */

  PetscInt cstart;  /* Contingency list start index for this processor */
  PetscInt cend;    /* Contingency list end idx (cstart+nc) for this processor */
		       
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

  PetscInt mode; /* Preventive or corrective mode */

  PetscBool iscoupling; /* Is each scenario coupled with base scenario? */
  PetscBool replicate_basecase; /* Replicate base case for all scenarios */

  void *solver; /* Solver object */
  struct _p_SCOPFLOWSolverOps solverops;
  char   solvername[64];

  /* List of solvers registered */
  struct _p_SCOPFLOWSolverList SCOPFLOWSolverList[SCOPFLOWSOLVERSMAX];
  PetscInt nsolversregistered;
  PetscBool SCOPFLOWSolverRegisterAllCalled;

  void    *solverdata;

  /* Data for contingencies */
  ContingencyList ctgclist;      /* List of contingencies */
  PetscBool       ctgcfileset;   /* Is the contingency file set ? */
  char            ctgcfile[100]; /* Contingency file */

  PetscBool       solutiontops;

  ContingencyFileInputFormat ctgcfileformat;
};

/* Register all SCOPFLOW solvers */
extern PetscErrorCode SCOPFLOWSolverRegisterAll(SCOPFLOW);
extern PetscErrorCode SCOPFLOWGetConstraints(SCOPFLOW,PetscInt,Vec*);
extern PetscErrorCode SCOPFLOWGetConstraintMultipliers(SCOPFLOW,PetscInt,Vec*);
extern PetscErrorCode SCOPFLOWReadContingencyData(SCOPFLOW,ContingencyFileInputFormat,const char[]);


#endif
