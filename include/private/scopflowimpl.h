
/**
 * @file scopflowimpl.h
 * @brief Private header file that defines data and structures for Security constrained Optimal Power Flow application
 */
#ifndef SCOPFLOWIMPL_H
#define SCOPFLOWIMPL_H

#include <ps.h>
#include <private/psimpl.h>
#include <opflow.h>
#include <tcopflow.h>
#include <scopflow.h>
#include <private/contingencylist.h>

#define SCOPFLOWSOLVERSMAX 10
#define SCOPFLOWMODELSMAX  10

struct _p_SCOPFLOWModelOps {
  PetscErrorCode (*destroy)(SCOPFLOW);
  PetscErrorCode (*setup)(SCOPFLOW);
  PetscErrorCode (*setnumvariablesandconstraints)(SCOPFLOW,PetscInt*,PetscInt*,PetscInt*,PetscInt*); /* Set number of variables and constraints for each contingency, and the number of equality/inequality coupling constraints */
  PetscErrorCode (*setvariablebounds)(SCOPFLOW,Vec,Vec); /* Upper and lower bounds on the vector */
  PetscErrorCode (*setconstraintbounds)(SCOPFLOW,Vec,Vec); /* Lower and upper bounds on constraints */
  PetscErrorCode (*setvariableandconstraintbounds)(SCOPFLOW,Vec,Vec,Vec,Vec); /* Lower and upper bounds on variables and constraints */
  PetscErrorCode (*setinitialguess)(SCOPFLOW,Vec); /* Set the initial guess for the optimization */
  PetscErrorCode (*computeconstraints)(SCOPFLOW,Vec,Vec);
  PetscErrorCode (*computejacobian)(SCOPFLOW,Vec,Mat);
  PetscErrorCode (*computehessian)(SCOPFLOW,Vec,Vec,Mat);
  PetscErrorCode (*computeobjandgradient)(SCOPFLOW,Vec,PetscScalar*,Vec); /* Objective and gradient routine */
  PetscErrorCode (*computeobjective)(SCOPFLOW,Vec,PetscScalar*); /* Objective */
  PetscErrorCode (*computegradient)(SCOPFLOW,Vec,Vec); /* Gradient of the objective function */
};

struct _p_SCOPFLOWModelList{
  char name[32]; /* Name of the model */
  PetscErrorCode (*create)(SCOPFLOW); /* Model creation routine */
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

struct _p_SCOPFLOWSolverList{
  char name[32]; /* Name of the solver */
  PetscErrorCode (*create)(SCOPFLOW); /* Solver object creation routine */
};

/**
 * @brief private struct for security optimal power flow application
*/
struct _p_SCOPFLOW{
  /* Sizes */
  PetscInt nc,Nc; /* Number of local and global (total) contingencies */
  PetscInt nx,Nx; /* Local and global (total) number of variables */
  PetscInt ncon,Ncon; /* Number of constraints */
  PetscInt nconeq,Nconeq; /* Local and global number of equality constraints */
  PetscInt nconineq, Nconineq; /* Local and global number of inequality constraints */
  PetscInt *nconeqcoup,*nconineqcoup;     /* Number of equality and inequality coupling constraints for each contingency */
  PetscInt nconcoup,Nconcoup; /* Number of coupling constraints */
  PetscInt *nxi; /* Number of variables for each contingency */
  PetscInt *ngi; /* Number of constraints for each contingency (includes coupling constraints) */
  PetscInt *xstarti; /* Starting location for the variables for contingency i in the big X vector */
  PetscInt *gstarti; /* Starting location for the constraints for contingency i in the big G vector */

  PetscInt cstart;  /* Contingency list start index for this processor */
  PetscInt cend;    /* Contingency list end idx (cstart+nc) for this processor */
		       
  COMM comm; /**< Communicator context */

  OPFLOW opflow0;   /* Base-case OPFLOW..maintained on each process to ease the calculations */

  OPFLOW *opflows; /* Array of optimal power flow application objects.
  		      Each processor creates nc objects, one for each 
		      contingency */

  TCOPFLOW *tcopflows; /* Array of multi-period optimal power flow application objects */

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

  Mat *Jcoup;  /* Jacobian for the coupling constraints (one per contingency) */
  Mat *JcoupT; /* Transpose of the coupling Jacobian (one per contingency)*/

  PetscBool iscoupling; /* Is each contingency coupled with base scenario? */

  PetscInt mode; /* 0 - preventive: all generation at PV and PQ buses is fixed. Only ref. bus(es) are allowed
		        to deviate from their base-case set points
		    1 - corrective: all generators are allowed to deviate from their base-case set points
		        limited by their 30. min ramp rates */

  void* model; /* Model object */
  struct _p_SCOPFLOWModelOps modelops;
  char modelname[64];

  struct _p_SCOPFLOWModelList SCOPFLOWModelList[SCOPFLOWMODELSMAX];
  PetscInt nmodelsregistered;
  PetscBool SCOPFLOWModelRegisterAllCalled;

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

  ContingencyFileInputFormat ctgcfileformat;

  PetscBool ismultiperiod; /* Are we doing a multi-period OPF? */

};

extern PetscErrorCode SCOPFLOWModelRegisterAll(SCOPFLOW);
extern PetscErrorCode SCOPFLOWSolverRegisterAll(SCOPFLOW);
extern PetscErrorCode SCOPFLOWGetConstraints(SCOPFLOW,PetscInt,Vec*);
extern PetscErrorCode SCOPFLOWGetConstraintMultipliers(SCOPFLOW,PetscInt,Vec*);
extern PetscErrorCode SCOPFLOWReadContingencyData(SCOPFLOW,ContingencyFileInputFormat,const char[]);
extern PetscErrorCode SCOPFLOWGetNumVariablesandConstraints(SCOPFLOW,PetscInt*,PetscInt*);

#endif
