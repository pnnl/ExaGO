/**
 * @file tcopflowimpl.h
 * @brief Private header file that defines data and structures for multi-period optimal Power Flow application
 */
#ifndef TCOPFLOWIMPL_H
#define TCOPFLOWIMPL_H

#include <ps.h>
#include <private/psimpl.h>
#include <opflow.h>
#include <tcopflow.h>
#include <private/contingencylist.h>

#define TCOPFLOWMODELSMAX  10
#define TCOPFLOWSOLVERSMAX 3


struct _p_TCOPFLOWModelOps {
  PetscErrorCode (*destroy)(TCOPFLOW);
  PetscErrorCode (*setup)(TCOPFLOW);
  PetscErrorCode (*setnumvariablesandconstraints)(TCOPFLOW,PetscInt*,PetscInt*,PetscInt*); /* Set number of variables for buses and branches, and total number of variables */
  PetscErrorCode (*setvariablebounds)(TCOPFLOW,Vec,Vec); /* Upper and lower bounds on the vector */
  PetscErrorCode (*setconstraintbounds)(TCOPFLOW,Vec,Vec); /* Lower and upper bounds on constraints */
  PetscErrorCode (*setvariableandconstraintbounds)(TCOPFLOW,Vec,Vec,Vec,Vec); /* Lower and upper bounds on variables and constraints */
  PetscErrorCode (*setinitialguess)(TCOPFLOW,Vec); /* Set the initial guess for the optimization */
  PetscErrorCode (*computeconstraints)(TCOPFLOW,Vec,Vec);
  PetscErrorCode (*computejacobian)(TCOPFLOW,Vec,Mat);
  PetscErrorCode (*computehessian)(TCOPFLOW,Vec,Vec,Mat);
  PetscErrorCode (*computeobjandgradient)(TCOPFLOW,Vec,PetscScalar*,Vec); /* Objective and gradient routine */
  PetscErrorCode (*computeobjective)(TCOPFLOW,Vec,PetscScalar*); /* Objective */
  PetscErrorCode (*computegradient)(TCOPFLOW,Vec,Vec); /* Gradient of the objective function */
};

struct _p_TCOPFLOWModelList{
  char name[32]; /* Name of the model */
  PetscErrorCode (*create)(TCOPFLOW); /* Model creation routine */
};

struct _p_TCOPFLOWSolverOps {
  PetscErrorCode (*destroy)(TCOPFLOW);
  PetscErrorCode (*setup)(TCOPFLOW);
  PetscErrorCode (*solve)(TCOPFLOW);
  PetscErrorCode (*getobjective)(TCOPFLOW,PetscReal*);
  PetscErrorCode (*getsolution)(TCOPFLOW,PetscInt,Vec*);
  PetscErrorCode (*getconvergencestatus)(TCOPFLOW,PetscBool*);
  PetscErrorCode (*getconstraints)(TCOPFLOW,PetscInt,Vec*);
  PetscErrorCode (*getconstraintmultipliers)(TCOPFLOW,PetscInt,Vec*);
};

struct _p_TCOPFLOWSolverList{
  char name[32]; /* Name of the solver */
  PetscErrorCode (*create)(TCOPFLOW); /* Solver object creation routine */
};

/**
 * @brief private struct for multi-period optimal power flow application
*/
struct _p_TCOPFLOW{
  /* Sizes */
  PetscInt nt,Nt; /* Number of local and global time-steps */
  PetscInt nx,Nx; /* Local and global (total) number of variables */
  PetscInt Ncon,ncon; /* Number of constraints */
  PetscInt Nconeq,nconeq; /* Local and global number of equality constraints */
  PetscInt nconineq, Nconineq; /* Local and global number of inequality constraints */
  PetscInt *nconineqcoup;      /* Number of inequality coupling constraints */
  PetscInt Nconcoup;  /* Number of coupling constraints between time-steps */
  PetscInt *nxi; /* Number of variables for each scenario */
  PetscInt *ngi; /* Number of constraints for each scenario (includes coupling constraints) */
  PetscInt *xstarti; /* Starting location for the variables for scenario i in the big X vector */
  PetscInt *gstarti; /* Starting location for the constraints for scenario i in the big G vector */

  
  COMM comm; /**< Communicator context */
  OPFLOW *opflows; /* Array of optimal power flow application objects.
  		      Each processor creates nt objects, one for each 
		      scenario */

  PetscBool setupcalled; /* TCOPFLOWSetUp called? */

  PetscBool converged; /* Convergence status */
  PetscInt  numiter; /* Number of iterations */
  PetscReal tolerance; /* Tolerance for TCOPFLOW */

  char netfile[PETSC_MAX_PATH_LEN]; /* Network data file */

  Vec  X,localX;    /* Global and local solution vector */
  Vec  G; /**< Inequality and equality constraint function */
  Vec  Lambda;

  Mat  Jac; /* Jacobian for equality and inequality constraints */

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

  void* model; /* Model object */
  struct _p_TCOPFLOWModelOps modelops;
  char modelname[64];

  struct _p_TCOPFLOWModelList TCOPFLOWModelList[TCOPFLOWMODELSMAX];
  PetscInt nmodelsregistered;
  PetscBool TCOPFLOWModelRegisterAllCalled;

  void *solver; /* Solver object */
  struct _p_TCOPFLOWSolverOps solverops;
  char   solvername[64];

  /* List of solvers registered */
  struct _p_TCOPFLOWSolverList TCOPFLOWSolverList[TCOPFLOWSOLVERSMAX];
  PetscInt nsolversregistered;
  PetscBool TCOPFLOWSolverRegisterAllCalled;

  void    *solverdata;

  PetscReal  dT; /* Time-step (in minutes) */
  PetscReal  duration; /* Time horizon (in hours) */

  /* Data for time-periods */
  char ploadprofile[PETSC_MAX_PATH_LEN]; /* Active load profile */
  PetscBool ploadprofileset; /* Is the active load power profile set? */
  char qloadprofile[PETSC_MAX_PATH_LEN]; /* Reactive load profile */
  PetscBool qloadprofileset; /* Is the reactive load power profile set ? */
  char windgenprofile[PETSC_MAX_PATH_LEN]; /* Wind generation profiles */
  PetscBool windgenprofileset; /* Is the wind generation profile set ? */

  Contingency *ctgc;
};

/* Register all TCOPFLOW solvers */
extern PetscErrorCode TCOPFLOWSolverRegisterAll(TCOPFLOW);
/* Register all TCOPFLOW models */
extern PetscErrorCode TCOPFLOWModelRegisterAll(TCOPFLOW);
extern PetscErrorCode TCOPFLOWReadPloadProfile(TCOPFLOW,char[]);
extern PetscErrorCode TCOPFLOWReadQloadProfile(TCOPFLOW,char[]);
extern PetscErrorCode TCOPFLOWReadWindGenProfile(TCOPFLOW,char[]);
extern PetscErrorCode TCOPFLOWGetConstraints(TCOPFLOW,PetscInt,Vec*);
extern PetscErrorCode TCOPFLOWGetConstraintMultipliers(TCOPFLOW,PetscInt,Vec*);
extern PetscErrorCode TCOPFLOWSetContingency(TCOPFLOW,Contingency*);


#endif
