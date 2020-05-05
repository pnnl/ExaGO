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

#define TCOPFLOWSOLVERSMAX 3

struct _p_TCOPFLOWSolverList{
  char name[32]; /* Name of the solver */
  PetscErrorCode (*create)(TCOPFLOW); /* Solver object creation routine */
};

struct _p_TCOPFLOWSolverOps {
  PetscErrorCode (*destroy)(TCOPFLOW);
  PetscErrorCode (*setup)(TCOPFLOW);
  PetscErrorCode (*solve)(TCOPFLOW);
};


/**
 * @brief private struct for multi-period optimal power flow application
*/
struct _p_TCOPFLOW{
  /* Sizes */
  PetscInt ns,Ns; /* Number of local and global time-steps */
  PetscInt nx,Nx; /* Local and global (total) number of variables */
  PetscInt Ncon,ncon; /* Number of constraints */
  PetscInt Nconeq,nconeq; /* Local and global number of equality constraints */
  PetscInt nconineq, Nconineq; /* Local and global number of inequality constraints */
  PetscInt *nconineqcoup;      /* Number of inequality coupling constraints */
  
  COMM comm; /**< Communicator context */
  OPFLOW *opflows; /* Array of optimal power flow application objects.
  		      Each processor creates ns objects, one for each 
		      scenario */

  PetscBool setupcalled; /* TCOPFLOWSetUp called? */

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

  void *solver; /* Solver object */
  struct _p_TCOPFLOWSolverOps solverops;

  /* List of solvers registered */
  struct _p_TCOPFLOWSolverList TCOPFLOWSolverList[TCOPFLOWSOLVERSMAX];
  PetscInt nsolversregistered;
  PetscBool TCOPFLOWSolverRegisterAllCalled;

  void    *solverdata;

  PetscReal  dT; /* Time-step (in minutes) */
  PetscReal  duration; /* Time horizon (in hours) */
  PetscInt   ntimesteps; /* Number of time-steps - This needs to replace Ns */

  /* Data for time-periods */
  char ploadprofile[100]; /* Active load profile */
  PetscBool ploadprofileset; /* Is the active load power profile set? */
  char qloadprofile[100]; /* Reactive load profile */
  PetscBool qloadprofileset; /* Is the reactive load power profile set ? */
  char windgenprofile[100]; /* Wind generation profiles */
  PetscBool windgenprofileset; /* Is the wind generation profile set ? */
  
};

/* Register all TCOPFLOW solvers */
extern PetscErrorCode TCOPFLOWSolverRegisterAll(TCOPFLOW);
extern PetscErrorCode TCOPFLOWReadPloadProfile(TCOPFLOW,char[]);
extern PetscErrorCode TCOPFLOWReadQloadProfile(TCOPFLOW,char[]);
extern PetscErrorCode TCOPFLOWReadWindGenProfile(TCOPFLOW,char[]);

#endif
