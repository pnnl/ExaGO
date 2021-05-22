
/**
 * @file sopflowimpl.h
 * @brief Private header file that defines data and structures for Stochastic Optimal Power Flow application
 */
#ifndef SOPFLOWIMPL_H
#define SOPFLOWIMPL_H

#include <ps.h>
#include <private/psimpl.h>
#include <opflow.h>
#include <sopflow.h>
#include <scopflow.h>
#include <private/scenariolist.h>

#define SOPFLOWSOLVERSMAX 10
#define SOPFLOWMODELSMAX  10

struct _p_SOPFLOWModelOps {
  PetscErrorCode (*destroy)(SOPFLOW);
  PetscErrorCode (*setup)(SOPFLOW);
  PetscErrorCode (*setnumvariablesandconstraints)(SOPFLOW,PetscInt*,PetscInt*,PetscInt*,PetscInt*); /* Set number of variables and the equality/inequality coupling constraints for each scenario */
  PetscErrorCode (*setvariablebounds)(SOPFLOW,Vec,Vec); /* Upper and lower bounds on the vector */
  PetscErrorCode (*setconstraintbounds)(SOPFLOW,Vec,Vec); /* Lower and upper bounds on constraints */
  PetscErrorCode (*setvariableandconstraintbounds)(SOPFLOW,Vec,Vec,Vec,Vec); /* Lower and upper bounds on variables and constraints */
  PetscErrorCode (*setinitialguess)(SOPFLOW,Vec); /* Set the initial guess for the optimization */
  PetscErrorCode (*computeconstraints)(SOPFLOW,Vec,Vec);
  PetscErrorCode (*computejacobian)(SOPFLOW,Vec,Mat);
  PetscErrorCode (*computehessian)(SOPFLOW,Vec,Vec,Mat);
  PetscErrorCode (*computeobjandgradient)(SOPFLOW,Vec,PetscScalar*,Vec); /* Objective and gradient routine */
  PetscErrorCode (*computeobjective)(SOPFLOW,Vec,PetscScalar*); /* Objective */
  PetscErrorCode (*computegradient)(SOPFLOW,Vec,Vec); /* Gradient of the objective function */
};

struct _p_SOPFLOWModelList{
  char name[32]; /* Name of the model */
  PetscErrorCode (*create)(SOPFLOW); /* Model creation routine */
};

struct _p_SOPFLOWSolverOps {
  PetscErrorCode (*destroy)(SOPFLOW);
  PetscErrorCode (*setup)(SOPFLOW);
  PetscErrorCode (*solve)(SOPFLOW);
  PetscErrorCode (*getobjective)(SOPFLOW,PetscReal*);
  PetscErrorCode (*getsolution)(SOPFLOW,PetscInt,Vec*);
  PetscErrorCode (*getconvergencestatus)(SOPFLOW,PetscBool*);
  PetscErrorCode (*getconstraints)(SOPFLOW,PetscInt,Vec*);
  PetscErrorCode (*getconstraintmultipliers)(SOPFLOW,PetscInt,Vec*);
};

struct _p_SOPFLOWSolverList{
  char name[32]; /* Name of the solver */
  PetscErrorCode (*create)(SOPFLOW); /* Solver object creation routine */
};

/**
 * @brief private struct for security optimal power flow application
*/
struct _p_SOPFLOW{
  /* Sizes */
  PetscInt ns,Ns; /* Number of local and global (total) scenarios */
  PetscInt nx,Nx; /* Local and global (total) number of variables */
  PetscInt ncon,Ncon; /* Number of constraints */
  PetscInt nconeq,Nconeq; /* Local and global number of equality constraints */
  PetscInt nconineq, Nconineq; /* Local and global number of inequality constraints */
  PetscInt *nconeqcoup, *nconineqcoup;     /* Number of equality/inequality coupling constraints for each scenario */
  PetscInt nconcoup,Nconcoup; /* Number of coupling constraints */
  PetscInt *nxi; /* Number of variables for each scenario */
  PetscInt *ngi; /* Number of constraints for each scenario (includes coupling constraints) */
  PetscInt *xstarti; /* Starting location for the variables for scenario i in the big X vector */
  PetscInt *gstarti; /* Starting location for the constraints for scenario i in the big G vector */

  PetscInt sstart;  /* Scenario list start index for this processor */
  PetscInt send;    /* Scenario list end idx (cstart+nc) for this processor */
		       
  COMM comm; /**< Communicator context */
  OPFLOW *opflows; /* Array of optimal power flow application objects.
  		      Each processor creates ns objects, one for each 
		      scenario */

  SCOPFLOW *scopflows; /* Array of security-constrained optimal power flow application objects. */
  PetscBool ismulticontingency; /* Is it a multi-contingency SOPFLOW? */


  PetscBool setupcalled; /* SOPFLOWSetUp called? */

  PetscBool converged; /* Convergence status */
  PetscInt  numiter; /* Number of iterations */
  PetscReal tolerance; /* Tolerance for SOPFLOW */

  char netfile[PETSC_MAX_PATH_LEN]; /* Network data file */

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

  PetscInt mode; /* 0 - preventive: all non-renewable generation at PV and PQ buses is fixed. Only generation at ref. bus(es) and renewable generation (wind, solar) are allowed to deviate from their base-case set points
		    1 - corrective: all generators are allowed to deviate from their base-case set points
		        limited by their 30. min ramp rates */

  void* model; /* Model object */
  struct _p_SOPFLOWModelOps modelops;
  char modelname[64];

  struct _p_SOPFLOWModelList SOPFLOWModelList[SOPFLOWMODELSMAX];
  PetscInt nmodelsregistered;
  PetscBool SOPFLOWModelRegisterAllCalled;

  void *solver; /* Solver object */
  struct _p_SOPFLOWSolverOps solverops;
  char   solvername[64];

  /* List of solvers registered */
  struct _p_SOPFLOWSolverList SOPFLOWSolverList[SOPFLOWSOLVERSMAX];
  PetscInt nsolversregistered;
  PetscBool SOPFLOWSolverRegisterAllCalled;

  void    *solverdata;

  /* Data for scenarios */
  ScenarioList    scenlist;
  PetscBool       scenfileset;   /* Is the scenario file set ? */
  char            scenfile[PETSC_MAX_PATH_LEN]; /* Scenario file */
  ScenarioFileInputFormat scenfileformat;
  ScenarioUncertaintyType scenunctype;

  MPI_Comm subcomm; /* Sub-communicators on which SCOPFLOW run */
};

extern PetscErrorCode SOPFLOWModelRegisterAll(SOPFLOW);
extern PetscErrorCode SOPFLOWSolverRegisterAll(SOPFLOW);
extern PetscErrorCode SOPFLOWGetConstraints(SOPFLOW,PetscInt,Vec*);
extern PetscErrorCode SOPFLOWGetConstraintMultipliers(SOPFLOW,PetscInt,Vec*);
extern PetscErrorCode SOPFLOWReadScenarioData(SOPFLOW,ScenarioFileInputFormat,const char[]);


#endif
