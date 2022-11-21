
/**
 * @file opflowimpl.h
 * @brief Private header file that defines data and structures for Optimal Power
 * Flow application
 */
#ifndef OPFLOWIMPL_H
#define OPFLOWIMPL_H

#include <opflow.h>
#include <pflow.h>
#include <private/psimpl.h>
#include <ps.h>

#define OPFLOWMODELSMAX 10
#define OPFLOWSOLVERSMAX 10

extern const char *const OPFLOWInitializationTypes[];

extern const char *const OPFLOWObjectiveTypes[];

extern const char *const OPFLOWGenBusVoltageTypes[];

typedef enum { DEFAULT = 0, HOST = 1, UM = 2, DEVICE = 3 } HIOPMemSpace;

extern const char *HIOPMemSpaceChoices[];

struct _p_OPFLOWModelOps {
  PetscErrorCode (*destroy)(OPFLOW);
  PetscErrorCode (*setup)(OPFLOW);
  PetscErrorCode (*setnumvariables)(
      OPFLOW, PetscInt *, PetscInt *,
      PetscInt *); /* Set number of variables for buses and branches, and total
                      number of variables */
  PetscErrorCode (*setnumconstraints)(
      OPFLOW, PetscInt *, PetscInt *, PetscInt *,
      PetscInt *); /* Set number of equality and inequality constraints */
  PetscErrorCode (*setvariablebounds)(
      OPFLOW, Vec, Vec); /* Upper and lower bounds on the vector */
  PetscErrorCode (*updatevariablebounds)(
      OPFLOW, Vec, Vec, void *); /* Upper and lower bounds on the vector */
  PetscErrorCode (*setvariableboundsarray)(
      OPFLOW, double *, double *); /* Array version of set variable bounds */
  PetscErrorCode (*setconstraintbounds)(
      OPFLOW, Vec, Vec); /* Lower and upper bounds on constraints */
  PetscErrorCode (*setconstraintboundsarray)(
      OPFLOW, double *, double *); /* Array version of constraint bounds */
  PetscErrorCode (*setvariableandconstraintbounds)(
      OPFLOW, Vec, Vec, Vec,
      Vec); /* Lower and upper bounds on variables and constraints */
  PetscErrorCode (*setvariableandconstraintboundsarray)(
      OPFLOW, double *, double *, double *,
      double *); /* Array version of setvariable and constraint bounds */
  PetscErrorCode (*setinitialguess)(
      OPFLOW, Vec, Vec); /* Set the initial guess for the optimization */
  PetscErrorCode (*setinitialguessarray)(
      OPFLOW, double *); /* Array version of set initial guess */
  PetscErrorCode (*computeequalityconstraints)(
      OPFLOW, Vec, Vec); /* Set equality constraints */
  PetscErrorCode (*computeequalityconstraintsarray)(
      OPFLOW, const double *,
      double *); /* Array version of set equality constraints */
  PetscErrorCode (*computeinequalityconstraints)(
      OPFLOW, Vec, Vec); /* Set inequality constraints */
  PetscErrorCode (*computeinequalityconstraintsarray)(
      OPFLOW, const double *, double *); /* Set equality constraints */
  PetscErrorCode (*computeconstraints)(OPFLOW, Vec, Vec);
  PetscErrorCode (*computeconstraintsarray)(
      OPFLOW, double *, double *); /* Array version of compute constraints */
  PetscErrorCode (*computeequalityconstraintjacobian)(OPFLOW, Vec, Mat);
  PetscErrorCode (*computeinequalityconstraintjacobian)(OPFLOW, Vec, Mat);
  PetscErrorCode (*computehessian)(OPFLOW, Vec, Vec, Vec, Mat);
  PetscErrorCode (*computeobjandgradient)(
      OPFLOW, Vec, PetscScalar *, Vec); /* Objective and gradient routine */
  PetscErrorCode (*computeobjective)(OPFLOW, Vec,
                                     PetscScalar *); /* Objective */
  PetscErrorCode (*computeobjectivearray)(
      OPFLOW, const double *,
      double *); /* Array version of compute objective */
  PetscErrorCode (*computegradient)(
      OPFLOW, Vec, Vec); /* Gradient of the objective function */
  PetscErrorCode (*computegradientarray)(
      OPFLOW, const double *, double *); /* Array version of compute gradient */
  PetscErrorCode (*computejacobian)(OPFLOW, Vec,
                                    Mat); /* Jacobian of the constraints */
  PetscErrorCode (*solutiontops)(
      OPFLOW); /* Update PS struct from OPFLOW solution */
  /* Following methods are only used with HIOP */
  PetscErrorCode (*computesparseequalityconstraintjacobianhiop)(
      OPFLOW, const double *, int *, int *,
      double *); /* Sparse Equality Constraints Jacobian */
  PetscErrorCode (*computesparseinequalityconstraintjacobianhiop)(
      OPFLOW, const double *, int *, int *,
      double *); /* Sparse Inequality Constraints Jacobian */

  PetscErrorCode (*computesparsehessianhiop)(OPFLOW, const double *,
                                             const double *, int *, int *,
                                             double *); /* Sparse Hessian */
  PetscErrorCode (*computedenseequalityconstraintjacobianhiop)(
      OPFLOW, const double *, double *); /* Dense Jacobian */
  PetscErrorCode (*computedenseinequalityconstraintjacobianhiop)(
      OPFLOW, const double *, double *); /* Dense Jacobian */
  PetscErrorCode (*computedensehessianhiop)(OPFLOW, const double *,
                                            const double *,
                                            double *); /* Dense Hessian */
  PetscErrorCode (*solutioncallbackhiop)(
      OPFLOW, const double *, const double *, const double *, const double *,
      const double *,
      double); // Call back for final solution

  /* Auxillary objective,gradient and hessian functions */
  /* Some applications may require to add custom objective function values in
     addition to what is being availble in OPFLOW. These auxillary functions
     provide setting custom objective, gradient, and hessian
  */
  PetscErrorCode (*computeauxobjective)(
      OPFLOW, const double *, double *,
      void *); /* Auxillary objective function */
  PetscErrorCode (*computeauxgradient)(OPFLOW, const double *, double *,
                                       void *); /* Auxillary gradient */
  PetscErrorCode (*computeauxhessian)(OPFLOW, const double *, Mat,
                                      void *); /* Auxillary hessian */
};

struct _p_OPFLOWSolverOps {
  PetscErrorCode (*destroy)(OPFLOW);
  PetscErrorCode (*setup)(OPFLOW);
  PetscErrorCode (*solve)(OPFLOW);
  PetscErrorCode (*getobjective)(OPFLOW, PetscReal *);
  PetscErrorCode (*getsolution)(OPFLOW, Vec *);
  PetscErrorCode (*getconvergencestatus)(OPFLOW, PetscBool *);
  PetscErrorCode (*getconstraints)(OPFLOW, Vec *);
  PetscErrorCode (*getconstraintmultipliers)(OPFLOW, Vec *);
};

struct _p_OPFLOWModelList {
  char name[32];                    /* Name of the model */
  PetscErrorCode (*create)(OPFLOW); /* Model creation routine */
};

struct _p_OPFLOWSolverList {
  char name[32];                    /* Name of the solver */
  PetscErrorCode (*create)(OPFLOW); /* Solver object creation routine */
};

/**
 * @brief private struct for optimal power flow application
 */
struct _p_OPFLOW {
  /* Common fields */
  COMM comm; /**< Communicator context */
  PS ps;     /**< Power system context */

  Vec X, localX; /* Global and local solution vector */
  Vec G;         /**< Inequality and equality constraint function */
  Vec Ge,
      Gelocal; /** < Equality constraint function vector (global and local) */
  Vec Gi;      /** < Inequality constraint function vector */

  Mat Jac;    /* Jacobian for equality and inequality constraints */
  Mat Jac_Ge; /* Equality constraint Jacobian */
  Mat Jac_Gi; /* Inequality constraint Jacobian */

  Mat Hes; /* Lagrangian Hessian */

  Vec Xl; /**< Lower bound on solution */
  Vec Xu; /**< Upper bound on solution */

  Vec Gl; /**< Lower bound on G */
  Vec Gu; /**< Upper bound on G */

  PetscScalar obj; /**< Objective function */
  Vec gradobj;     /**< Gradient of the objective function */

  PetscScalar obj_factor; /* IPOPT scales the objective hessian part with this
                             factor. For all other solvers, unless it is set,
                             obj_factor = 1.0. */

  PetscBool obj_gencost; /* Objective is generator cost minimization */

  PetscBool has_gensetpoint; /* Has a set-point for the real power generation */

  PetscBool use_agc; /* Uses AGC for generator dispatch */

  PetscBool setupcalled; /* OPFLOWSetUp called? */

  PetscReal tolerance; /* Tolerance for OPFLOW */

  OPFLOWInitializationType initializationtype; /* OPFLOW Initialization type */

  OPFLOWObjectiveType objectivetype; /* OPFLOW objective */

  OPFLOWGenBusVoltageType genbusvoltagetype; /* OPFLOW Genbus voltage type */

  PetscInt nconeq, Nconeq; /* Local and global number of equality constraints,
                              excluding ghosts! */
  PetscInt nconineq,
      Nconineq;        /* Local and global number of inequality constraints */
  PetscInt Ncon, ncon; /* Total number of constraints (equality + inequality) */
  PetscInt nx,
      Nx; /* Total number of local and global variables, excluding ghosts! */

  PetscInt nxsparse, nxdense; /* Only used by HIOP MDS */

  //  Mat Jac_GeT; /* Transpose of equality constraint Jacobian */
  //  Mat Jac_GiT; /* Transpose of inequality constraint Jacobian */

  /* Lagrange multipliers */
  Vec Lambda, Lambdae, Lambdai; // Lagrange multipliers for constraints Lambda =
                                // [Lambdae; Lambdai];
  Vec Lambdaelocal;             /* Local Lambdae vector */

  void *solver; /* Solver object */
  struct _p_OPFLOWSolverOps solverops;
  char solvername[64];

  void *model; /* Model object */
  struct _p_OPFLOWModelOps modelops;
  char modelname[64];

  char _p_hiop_compute_mode[64];
  int _p_hiop_verbosity_level;

  /* List of models and solvers registered */
  struct _p_OPFLOWModelList OPFLOWModelList[OPFLOWMODELSMAX];
  struct _p_OPFLOWSolverList OPFLOWSolverList[OPFLOWSOLVERSMAX];
  PetscInt nmodelsregistered;
  PetscInt nsolversregistered;
  PetscBool OPFLOWModelRegisterAllCalled;
  PetscBool OPFLOWSolverRegisterAllCalled;

  PetscBool ignore_lineflow_constraints; /* Ignore line flow constraints */

  PetscBool include_loadloss_variables; /* Include variables for loadloss */
  PetscReal loadloss_penalty;

  PetscBool include_powerimbalance_variables; /* Include variables for power
                                                 imbalance */
  PetscReal powerimbalance_penalty;

  PetscBool has_powersetpoint; /* Use real power set-point */

  PetscInt numits; /* Number of solver iterations */

  PetscInt nlinekvmon;    /* Number of line kv levels to monitor */
  PetscScalar *linekvmon; /* Line kv levels to monitor */

  PetscInt nlinesmon; /* Number of monitored lines */
  PetscInt *linesmon; /* List of lines monitored */

  /* Flag to denote if the OPFLOW solution has been transfered to PS struct via
   * OPLOWSolutionToPS call */
  PetscBool solutiontops;

  /* Global indices for the equality constraints. It is used when operating on
   * equality constraints */
  PetscInt *eqconglobloc;

  VecScatter
      scattereqcon; /* Vector scatter object for scattering from global equality
                       constraints to local equality constraints vector */

  /* Local and global index sets used in VecScatter */
  IS isconeqlocal;
  IS isconeqglob;

  /* Used only when initialization from power flow (ACPF) is called */
  PFLOW initpflow; /**< power flow solver context for obtaining initial
                      conditions */
  PetscSection initpflowpsection; /** < PetscSection object to hold dofs at each
                                     vertex */
  PetscSection initpflowpglobsection; /** <Global section for initpflow */
  PetscSection
      defaultsection; /** < PetscSection used with opflow. This is temporarily
                         stored when running the initial power flow */
  PetscSection defaultglobalsection; /** <Global section used with opflow. This
                                        is temporarily stored when running the
                                        initial power flow */
  DM initpflowdm;                    /* DM used with initial power flow */

  PetscInt *busnvararray,
      *branchnvararray; /* Number of variables at buses and branches */

  PetscBool spdnordering; /* TRUE if model uses sparse dense ordering */
  PetscInt *idxn2sd_map;  /* Mapping from natural to sparse-dense ordering */

  /** @brief Logging events that apply to interface */
  PetscLogEvent objlogger, gradlogger, eqconslogger, ineqconslogger,
      eqconsjaclogger, ineqconsjaclogger, hesslogger, solvelogger,
      densehesslogger, sparsehesslogger, denseineqconsjaclogger,
      denseeqconsjaclogger;

  /** @brief number of nonzeros used with HiOP MDS */
  PetscInt nnz_eqjacsp, nnz_ineqjacsp, nnz_hesssp;

  /** @brief user provided data struct for auxillary objective */
  void *userctx;

  PetscBool skip_options; /* Skip run-time options */

  PetscScalar weight; /* Weight for this system condition (0,1) */

  HIOPMemSpace mem_space; /* Memory space used with HIOP */
};

/* Registers all the OPFLOW models */
extern PetscErrorCode OPFLOWModelRegisterAll(OPFLOW);

/* Register all OPFLOW solvers */
extern PetscErrorCode OPFLOWSolverRegisterAll(OPFLOW);

/* Internal function to check model + solver compatibility */
extern PetscErrorCode OPFLOWCheckModelSolverCompatibility(OPFLOW);

/* Internal function to set number of variables */
extern PetscErrorCode OPFLOWSetNumVariables(OPFLOW, PetscInt *, PetscInt *,
                                            PetscInt *);

/* Internal function to set number of constraints */
extern PetscErrorCode OPFLOWSetNumConstraints(OPFLOW, PetscInt *, PetscInt *,
                                              PetscInt *, PetscInt *);

extern PetscErrorCode OPFLOWNaturalToSpDense(OPFLOW, const double *, double *);

extern PetscErrorCode OPFLOWSpDenseToNatural(OPFLOW, const double *, double *);

extern PetscErrorCode OPFLOWSetHIOPMemSpace(OPFLOW, HIOPMemSpace);
extern PetscErrorCode OPFLOWGetHIOPMemSpace(OPFLOW, HIOPMemSpace *);

#endif
