
/**
 * @file sopflowimpl.h
 * @brief Private header file that defines data and structures for Stochastic
 * Optimal Power Flow application
 */
#ifndef SOPFLOWIMPL_H
#define SOPFLOWIMPL_H

#include <opflow.h>
#include <private/psimpl.h>
#include <private/scenariolist.h>
#include <private/opflowimpl.h>
#include <ps.h>
#include <scopflow.h>
#include <sopflow.h>

#define SOPFLOWSOLVERSMAX 10
#define SOPFLOWMODELSMAX 10

struct _p_SOPFLOWModelOps {
  PetscErrorCode (*destroy)(SOPFLOW);
  PetscErrorCode (*setup)(SOPFLOW);
  PetscErrorCode (*setnumvariablesandconstraints)(
      SOPFLOW, PetscInt *, PetscInt *, PetscInt *,
      PetscInt *); /* Set number of variables and the equality/inequality
                      coupling constraints for each scenario */
  PetscErrorCode (*setvariablebounds)(
      SOPFLOW, Vec, Vec); /* Upper and lower bounds on the vector */
  PetscErrorCode (*setconstraintbounds)(
      SOPFLOW, Vec, Vec); /* Lower and upper bounds on constraints */
  PetscErrorCode (*setvariableandconstraintbounds)(
      SOPFLOW, Vec, Vec, Vec,
      Vec); /* Lower and upper bounds on variables and constraints */
  PetscErrorCode (*setinitialguess)(
      SOPFLOW, Vec); /* Set the initial guess for the optimization */
  PetscErrorCode (*computeconstraints)(SOPFLOW, Vec, Vec);
  PetscErrorCode (*computejacobian)(SOPFLOW, Vec, Mat);
  PetscErrorCode (*computehessian)(SOPFLOW, Vec, Vec, Mat);
  PetscErrorCode (*computeobjandgradient)(
      SOPFLOW, Vec, PetscScalar *, Vec); /* Objective and gradient routine */
  PetscErrorCode (*computebaseobjective)(
      SOPFLOW, Vec, PetscScalar *); /* Base case Objective */
  PetscErrorCode (*computetotalobjective)(SOPFLOW, Vec,
                                          PetscScalar *); /* Total Objective */

  PetscErrorCode (*computegradient)(
      SOPFLOW, Vec, Vec); /* Gradient of the objective function */
};

struct _p_SOPFLOWModelList {
  char name[32];                     /* Name of the model */
  PetscErrorCode (*create)(SOPFLOW); /* Model creation routine */
};

struct _p_SOPFLOWSolverOps {
  PetscErrorCode (*destroy)(SOPFLOW);
  PetscErrorCode (*setup)(SOPFLOW);
  PetscErrorCode (*solve)(SOPFLOW);
  PetscErrorCode (*getbaseobjective)(SOPFLOW, PetscReal *);
  PetscErrorCode (*gettotalobjective)(SOPFLOW, PetscReal *);
  PetscErrorCode (*getsolution)(SOPFLOW, PetscInt, Vec *);
  PetscErrorCode (*getconvergencestatus)(SOPFLOW, PetscBool *);
  PetscErrorCode (*getconstraints)(SOPFLOW, PetscInt, Vec *);
  PetscErrorCode (*getconstraintmultipliers)(SOPFLOW, PetscInt, Vec *);
};

struct _p_SOPFLOWSolverList {
  char name[32];                     /* Name of the solver */
  PetscErrorCode (*create)(SOPFLOW); /* Solver object creation routine */
};

/**
 * @brief private struct for security optimal power flow application
 */
struct _p_SOPFLOW {
  /* Sizes */
  PetscInt ns, Ns;         /* Number of local and global (total) scenarios */
  PetscInt nx, Nx;         /* Local and global (total) number of variables */
  PetscInt ncon, Ncon;     /* Number of constraints */
  PetscInt nconeq, Nconeq; /* Local and global number of equality constraints */
  PetscInt nconineq,
      Nconineq; /* Local and global number of inequality constraints */
  PetscInt *nconeqcoup, *nconineqcoup; /* Number of equality/inequality coupling
                                          constraints for each scenario */
  PetscInt nconcoup, Nconcoup;         /* Number of coupling constraints */
  PetscInt *nxi; /* Number of variables for each scenario */
  PetscInt *ngi; /* Number of constraints for each scenario (includes coupling
                    constraints) */
  PetscInt *xstarti; /* Starting location for the variables for scenario i in
                        the big X vector */
  PetscInt *gstarti; /* Starting location for the constraints for scenario i in
                        the big G vector */

  PetscInt sstart; /* Scenario list start index for this processor */
  PetscInt send;   /* Scenario list end idx (cstart+nc) for this processor */

  COMM comm;       /**< Communicator context */
  OPFLOW *opflows; /* Array of optimal power flow application objects.
                      Each processor creates ns objects, one for each
                      scenario */

  SCOPFLOW *scopflows; /* Array of security-constrained optimal power flow
                          application objects. */
  PetscBool ismulticontingency; /* Is it a multi-contingency SOPFLOW? */

  PetscBool ismultiperiod;               /* Is the contingency multi-period? */
  PetscInt Nc;                           /* Number of contingencies */
  PetscReal dT;                          /* Time-step */
  PetscReal duration;                    /* duration */
  char ploadprofile[PETSC_MAX_PATH_LEN]; /* Active load profile */
  char qloadprofile[PETSC_MAX_PATH_LEN]; /* Reactive load profile */

  OPFLOWInitializationType
      initialization_type; /* Initialization type for OPFLOW */
  OPFLOWGenBusVoltageType
      gen_bus_voltage_type; /* Gen bus voltage type for OPFLOW */
  PetscBool
      ignore_lineflow_constraints; /* Line flow constraints flag for OPFLOW */

  PetscBool flatten_contingencies; /* Flattens the contingencies to fuse the
                                      scenarios and contingencies */

  PetscBool setupcalled; /* SOPFLOWSetUp called? */

  PetscBool converged; /* Convergence status */
  PetscInt numiter;    /* Number of iterations */
  PetscReal tolerance; /* Tolerance for SOPFLOW */

  char netfile[PETSC_MAX_PATH_LEN]; /* Network data file */
  char ctgfile[PETSC_MAX_PATH_LEN]; /* Contingency data file */
  char windgen[PETSC_MAX_PATH_LEN]; /* Wingen data file */

  Vec X, localX; /* Global and local solution vector */
  Vec G;         /**< Inequality and equality constraint function */
  Vec Ge,
      Gelocal; /** < Equality constraint function vector (global and local) */
  Vec Gi;      /** < Inequality constraint function vector */
  Vec Lambdai;
  Vec Lambdae, Lambdaelocal;
  Vec Lambda;

  Mat Jac;    /* Jacobian for equality and inequality constraints */
  Mat Jac_Ge; /* Equality constraint Jacobian */
  Mat Jac_Gi; /* Inequality constraint Jacobian */

  Mat Hes; /* Lagrangian Hessian */

  Vec Xl; /**< Lower bound on solution */
  Vec Xu; /**< Upper bound on solution */

  Vec Gl; /**< Lower bound on G */
  Vec Gu; /**< Upper bound on G */

  PetscScalar objbase; /**< Base-case objective function value */
  PetscScalar objtot;  /**< Total objective function value */

  Vec gradobj; /**< Gradient of the objective function */

  PetscScalar obj_factor; /* IPOPT scales the objective hessian part with this
                             factor. For all other solvers, unless it is set,
                             obj_factor = 1.0. */

  Mat *Jcoup;  /* Jacobian for the coupling constraints (one per scenario) */
  Mat *JcoupT; /* Transpose of the coupling Jacobian (one per scenario)*/

  PetscBool iscoupling; /* Is each scenario coupled with base scenario? */

  PetscInt mode; /* 0 - preventive: all non-renewable generation at PV and PQ
                    buses is fixed. Only generation at ref. bus(es) and
                    renewable generation (wind, solar) are allowed to deviate
                    from their base-case set points 1 - corrective: all
                    generators are allowed to deviate from their base-case set
                    points limited by their 30. min ramp rates */

  void *model; /* Model object */
  struct _p_SOPFLOWModelOps modelops;
  char modelname[PETSC_MAX_PATH_LEN];

  struct _p_SOPFLOWModelList SOPFLOWModelList[SOPFLOWMODELSMAX];
  PetscInt nmodelsregistered;
  PetscBool SOPFLOWModelRegisterAllCalled;

  void *solver; /* Solver object */
  struct _p_SOPFLOWSolverOps solverops;
  char solvername[PETSC_MAX_PATH_LEN];

  /* List of solvers registered */
  struct _p_SOPFLOWSolverList SOPFLOWSolverList[SOPFLOWSOLVERSMAX];
  PetscInt nsolversregistered;
  PetscBool SOPFLOWSolverRegisterAllCalled;

  void *solverdata;

  /* Data for scenarios */
  ScenarioList scenlist;
  PetscBool scenfileset;             /* Is the scenario file set ? */
  char scenfile[PETSC_MAX_PATH_LEN]; /* Scenario file */
  ScenarioFileInputFormat scenfileformat;
  ScenarioUncertaintyType scenunctype;

  MPI_Comm subcomm; /* Sub-communicators on which SCOPFLOW run */

  OPFLOW opflow0; /* Base scenario, each rank has this information */

  /* Used when flatten contingency option is chosen */
  /* These are arrays of size ns on each rank where ns is the number of local
     cases A case is a scenario-contingency combination
  */
  PetscInt *scen_num; /* scenario number */
  PetscInt *cont_num; /* contingency number */
  /* Data for contingencies */
  ContingencyList ctgclist;          /* List of contingencies */
  PetscBool ctgcfileset;             /* Is the contingency file set ? */
  char ctgcfile[PETSC_MAX_PATH_LEN]; /* Contingency file */
  ContingencyFileInputFormat ctgcfileformat;

  // Only used for HIOP solver
  std::string subproblem_model;
  std::string subproblem_solver;
  std::string compute_mode;
  int verbosity_level;
  std::string mem_space; /* Memory space used with HIOP */

  /** @brief Logging events that apply to interface */
  PetscLogEvent outputlogger;
};

extern PetscErrorCode SOPFLOWModelRegisterAll(SOPFLOW);
extern PetscErrorCode SOPFLOWSolverRegisterAll(SOPFLOW);
extern PetscErrorCode SOPFLOWGetConstraints(SOPFLOW, PetscInt, Vec *);
extern PetscErrorCode SOPFLOWGetConstraintMultipliers(SOPFLOW, PetscInt, Vec *);
extern PetscErrorCode SOPFLOWReadScenarioData(SOPFLOW, ScenarioFileInputFormat,
                                              const char[]);

#endif
