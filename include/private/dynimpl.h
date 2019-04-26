/**
 * @file dynimpl.h
 * @brief Private header file for dynamics simulation
 */

#ifndef DYNIMPL_H
#define DYNIMPL_H

#include <ps.h>
#include <pflow.h>
#include <dyn.h>

#define MAXDYNEVENTTYPES 10 /**< Maximum number of distinct event types allowed */

#define MAXDYNEVENTS 100 /**< Max. allowed events */

#define MAXDYNEVENTSLOCATEDATSTEP 100 /**< Maximum number of event locations at a single step */

typedef struct _p_DYNEvent *DYNEvent;

extern PetscInt ndyneventtypesregistered; /**< Number of event types registered */

/**
 * @brief private struct to manage the list of registered events
 */
struct _p_DYNEventTypesList{
  char name[16]; /**< Name of the model (String identifier) */
  PetscErrorCode (*create)(DYNEvent); /**< Class constructor */
};

extern struct _p_DYNEventTypesList DYNEventTypesList[MAXDYNEVENTTYPES];

typedef const char* DYNEventType;

typedef enum{
  DYNEVENT_LOCATION_NOTSET,
  DYNEVENT_ON_BUS,
  DYNEVENT_ON_LINE,
  DYNEVENT_ON_EXC
}DYNEventLocation;

typedef enum{FREQVIOL,FREQ}DYNOBJFUNTYPE;

#define DYNEventFault "'FAULT'"
#define DYNEventLineSwitching "'LINESW'"
#define DYNEventGenTrip "'GENTRIP'"
#define DYNEventExciterLimits "'EXCLIM'"

/** 
 * @brief The event context
 */
struct _p_DYNEvent{
  char name[16];  /**< Name of the event (String identifier in the file) */
  void *data; /**< Private implementation data */
  PetscInt numMonitors; /**< Number of monitors */
  DYNEventLocation location; /**< Location of the event (bus or line) */
  PSBUS bus; /**< The bus structure */
  PSLINE line; /**< The line structure */
  PetscInt  *direction; /**< direction flags for event monitors */
  PetscBool *terminate; /**< terminate flags for event monitors */
  char *eventstring;
  /* Operations on DYNEvent */
  PetscErrorCode (*readdata)(DYNEvent,char*); /**< Read event data */
  PetscErrorCode (*getbusnumber)(DYNEvent,PetscInt*); /**< Get bus number */
  PetscErrorCode (*getfromtobusnumbers)(DYNEvent,PetscInt*,PetscInt*); /**< Get from and to bus numbers (for events on line) */
  PetscErrorCode (*eventmonitor)(DYNEvent,PetscReal,Vec,PetscScalar*); /**< Event monitoring function */
  PetscErrorCode (*posteventfunction)(DYNEvent,PetscInt,PetscInt[],PetscReal,Vec,PetscBool,PetscBool*); /**< Post event handling routine */
  PetscErrorCode (*posteventdgamma)(DYNEvent,PetscInt,PetscInt[],PetscReal,Vec,PetscBool); /**< Post event routine for sensitivity calculations */
  PetscErrorCode (*destroy)(DYNEvent); /**< Destroy the event object */
};

# define MAXPOSTSTEPFCNS 10
# define MAXPRESTEPFCNS  10  

/**
 * @brief private dynamic simulation application struct
 */
struct _p_DYN{
  COMM comm;   /**< Communicator context */
  PS   ps;     /**< Power system context */
  
  Vec X; /**< Solution vector */
  
  Mat Jac; /**< Jacobian matrix */
  
  TS ts; /**< The time-stepping solver */

  PFLOW initpflow; /**< power flow solver context for obtaining initial conditions */
  PetscSection initpflowpsection; /** < PetscSection object to hold dofs at each vertex */
  PetscSection initpflowpglobsection; /** <Global section for initpflow */
  PetscSection defaultsection;      /** < PetscSection used with dyn. This is temporarily stored when running the initial power flow */
  PetscSection defaultglobalsection; /** <Global section used with dyn. This is temporarily stored when running the initial power flow */
  DM           initpflowdm; /* DM used with initial power flow */

  PetscInt ndiff; /**< Number of differential equations */
  IS isdiff;    /**< Index set to hold only the differential variables */
  PetscBool solve_alg_only; /**< Flag to indicate that only the algebraic equations are being solved. This
                               is TRUE when an event (for example) causes discontinuity and the algebraic
                               equations need to be resolved to obtain consistent initial conditions*/
  PetscBool setupcalled; /**< DYNSetUp called? */
  PetscBool prepflow; /**< Flag to indicate if a power flow needs to be run before commencing dynamics simulation */
  PetscBool prepflowconverged; /**< Flag to indicate if a power flow converges */

  PetscInt npoststepfcns; /**<Number of user post-step callback functions */
  PetscErrorCode (*poststepfcn[MAXPOSTSTEPFCNS])(DYN);

  PetscInt nprestepfcns; /**<Number of user pre-step callback functions */
  PetscErrorCode (*prestepfcn[MAXPRESTEPFCNS])(DYN);

  /* Data for DYNEvent */
  PetscInt Nevents;           /**< Number of events on all processors */
  PetscInt nevents;           /**< Number of events on this processor */
  PetscInt toteventmonitors;    /**< Total number of monitors (event functions) for events */
  PetscInt event_idx[MAXDYNEVENTS]; /**< The locations of the events for this processor in the events array */
  DYNEvent events[MAXDYNEVENTS]; /**< Array of events on this processor */
  PetscInt *eventdirection;/**< directions for zero crossing..see TSSetEventMonitor */
  PetscBool *eventterminate; /**< termination flags..see TSSetEventMonitor */

  /* Time duration and step information */
  PetscInt  t0;       /**< Start time */
  PetscReal tmax;     /**< Max. time */
  PetscInt  maxsteps; /**< Max. integration steps */
  PetscReal dt0;      /**< First time step */

  Vec       vatol,vrtol; /**< Vector of tolerances for time-step adapativity */
  /* Activate for using semi-explicit scheme */
  PetscBool use_semiexplicit; /* Option -dyn_use_semiexplicit */

  /* For sensitivity calculation */
  Vec *lambda;  /**< Gradients of cost functions w.r.t initial conditions */
  Vec *mu;      /**< Gradients of cost functions w.r.t parameters */
  Vec  *dy0dp;   /**< Partial derivatives of initial conditions w.r.t parameters */
  Mat  Jacp;    /**< Jacobian of the DAE RHS w.r.t parameters */
  Mat  ICp;    /**< Jacobian of the Initial Conditions w.r.t parameters */
  Vec  costintegral; /**< Cost Integral */
#ifdef FWDSA
  Vec  *s,*sp;
  Vec  *intgrad,*intgradp;
  Vec  *fwdJacp;
#endif
  /* Cost function
     r = scal*f^exp -> where f is the cost function and r is the constraint
     that gets incorporated in the optimization
  */
  PetscReal   scal; /**< Scaling of the cost function */
  PetscReal   exp; /**< Exponent for smoothing cost function */
  PetscReal   eta; /**< Constraint bound (r <= \eta) */

  PetscScalar freq_min; /**< Min. frequency limits */
  PetscScalar freq_max; /**< Max. frequency limits */
  PetscBool   sum_freq_costfun;  /**< Sum all frequency cost functions */
  PetscInt    nparams; /**< Number of parameters (local) */
  PetscInt    Nparams; /**< Number of parameters (global) */
  PetscInt    ncostfcns; /**< Number of sensitivity cost functions */
  PetscScalar *costfcns; /**< Values for cost functions, memory allocated by users */
  Vec         Xparam;   /**< Temporary work vector (size = number of parameters) */

  PetscBool useadjoint;
  PetscBool useforward;
  PetscBool usefinitediff;
  PetscBool printsensi;
  PetscBool monitor;
  Vec       X0;
  Vec       Xtmp; /**< Temporary work vector duplicated from X*/
  Vec       Xtmp2; /**< Another temporary work vector for jump conditions */
  Vec       Xtmp3; /**< Another temporary work vector for costintegrals */
  Vec       dGammadX;
  Vec       dGammadP;
  PetscBool eval_integral;
  DYNOBJFUNTYPE   objfun;
  PetscScalar *feventtmp; /**< Temporary array used for checking and setting machine limits on disturbance on/off times */
  char   viewer_dir[PETSC_MAX_PATH_LEN]; /**< Directory location where the output files for visualization are saved */
  PetscBool visualize; /**< flag set if visualization front end is used */
  PetscLogEvent  logdynsolve; /**< Logging object for DYN time-stepping */
  PetscInt eventlogcount; /**< Event Log Counter (SAM)*/
  PetscBool tsconverged;
  PetscReal finaltime; /**< The actual time that TS solvers stops. TS may fail to converge in the middle. */

  void*     userdata;

  PetscBool switching_solution; /* Flag to indicate if a switching solution is being performed (only the solution of algebraic variables) 
				  Its set during algebraic solve and reset during prestep */
};

/**
 * @brief Registers all event types
 * Notes: This routine is called only once to set up DYNEventTypesList
 */
extern PetscErrorCode DYNEventRegisterAll(void);
/**
 * @brief Creates an instance of DYNEvent
 * @param [out] DYNEvent* dyneventout - the DYNEvent object
 */
extern PetscErrorCode DYNEventCreate(DYNEvent*);
/**
 * @brief Destroy an instance of DYNEvent 
 * @param [in] DYNEvent* dyneventout - the DYNEvent object
 */
extern PetscErrorCode DYNEventDestroy(DYNEvent*);
/**
 * @brief Sets the info for monitoring the event
 * @param [in] DYNEvent dyneventout - the DYNEvent object
 * @param [in] PetscInt numMonitors - number of monitors for this event
 * @param [in] PetscInt[] direction   - directions for zero crossing detection for each monitor ( 1 -> going positive, -1 -> going negative, 0 for both directions)
 * @param [in] PetscBool[] terminate   - termination flags for each monitor if one wishes to terminate the time-stepping after an event is detected
 * @param [in] PetscErrorCode (*eventmonitor)(DYNEvent,PetscReal,Vec,PetscScalar*) eventmonitor - the monitoring function for this event. This will be called after each time-step
 * @param [in] PetscErrorCode (*posteventfunction)(DYNEvent,PetscInt,PetscInt[],PetscReal,Vec,PetscBool,PetscBool*) posteventfunction - optional postevent function to allow taking actions after an event has occured. This routine is called after an event has occured
 */
extern PetscErrorCode DYNEventSetMonitors(DYNEvent,PetscInt,PetscInt[],PetscBool[],PetscErrorCode (*eventmonitor)(DYNEvent,PetscReal,Vec,PetscScalar*),PetscErrorCode (*posteventfunction)(DYNEvent,PetscInt,PetscInt[],PetscReal,Vec,PetscBool,PetscBool*));
/**
 * @brief NOT IMPLIMENTED
 */
extern PetscErrorCode DYNEventAddMonitors(DYNEvent,PetscInt,PetscInt[],PetscBool[]);
/**
 * @brief NO COMMENT
 */
extern PetscErrorCode DYNSetUpEvents(DYN);
/**
 * @brief NO COMMENT
 */
extern PetscErrorCode DYNEventSetLocation(DYNEvent,DYNEventLocation);
/**
 * @brief Called by TS after every step
 * @param [in] TS ts - the TS object
 * @param [in] PetscReal t  - the current time
 * @param [unknown] Vec X
 * @param [in] PetscScalar* f  - array of event function residuals
 * @param [in] void* ctx - application specific context
 * Notes: DYNEventMonitor() is also called by the post event function to set the feventtmp array. This array is used in checking if any machine device limits have been violated.
 */
extern PetscErrorCode DYNEventMonitor(TS,PetscReal,Vec,PetscScalar*,void*);
/**
 * @brief Post event function routine
 * @param [unknown] TS ts
 * @param [unknown] PetscInt nmonitors
 * @param [unknown] PetscInt[] monitor_list
 * @param [unknown] PetscReal t
 * @param [unknown] Vec X
 * @param [unknown] PetscBool forwardsolve
 * @param [unknown] void* ctx
 */
extern PetscErrorCode DYNEventPostFunction(TS,PetscInt,PetscInt[],PetscReal,Vec,PetscBool,void*);
/**
 * @brief Residual function for SNES to solve the algebraic part when an event is triggered
 * @param [unknown] SNES snes
 * @param [unknown] Vec X
 * @param [unknown] Vec F
 * @param [unknown] void* ctx
 */
extern PetscErrorCode DYNEventPostSolveAlgebraic(SNES,Vec,Vec,void*);
/**
 * @brief Adds shunt conductance and susceptance at the incident bus for this dynevent
 * @param [in] DYNEvent dynevent - The DYNEvent object
 * @param [in] PetscScalar Gs       - shunt conductance
 * @param [in] PetscScalar Bs       - shunt susceptance
 * Notes: Gs and Bs should be in per unit
 */
extern PetscErrorCode DYNEventAddBusShunt(DYNEvent,PetscScalar,PetscScalar);

/**
 * @brief Sets the generator status at a bus
 * @param [in] DYNEvent dynevent - The DYNEvent object
 * @param [in] PetscInt gbus       - Generator bus number
 * @param [in] char[]   gid        - Generator id
 * @param [in] int      status     - Generator status
 */
extern PetscErrorCode DYNEventSetGenStatus(DYNEvent,PetscInt,char[],PetscInt);

/**
 * @brief Sets the status of the line associated with this event
 * @param [in] DYNEvent dynevent - The DYNEvent object
 * @param [in] PetscInt status   - (0 or 1) the status of the line
 */
extern PetscErrorCode DYNEventSetLineStatus(DYNEvent,PetscInt);

/**
 * @brief Sets the events on machines and their control circuity
 * @param [in] DYN dyn - the DYN object
 * Notes: This routine adds events for gens, excs, and other controllers.
 */
extern PetscErrorCode DYNSetMachineEvents(DYN);
/**
 * @brief Sets the events on machines and their control circuity
 * @param [in] DYN dyn - the DYN object
 * @param [in] PetscReal t   - the current time
 * @param [in] Vec X   - the solution vector at current time
 * @param [unknown] PetscBool forwardsolve
 */

/**
 * @brief The IFunction for the time-stepping solver
 * @param [unknown] TS ts
 * @param [unknown] PetscReal t  
 * @param [unknown] Vec X 
 * @param [unknown] Vec Xdot 
 * @param [unknown] Vec F 
 * @param [unknown] void* ctx
 * Notes:
 * See PETSc routine TSSetIFunction()\n
 * The equations are expressed in current balance form:\n
 * I_gen - I_network (=Y*V) + I_shunt - I_load = 0\n
 * The current balance equations for each bus are ordered as [I_bus(imag);I_bus(real)]
 * Ordering them in such a manner allows having the susceptance B in the diagonal location
 * of the Jacobian matrix. For transmission networks, B > G (in many cases G=0) and hence
 * having this ordering is better.
 */
extern PetscErrorCode DYNIFunction(TS,PetscReal,Vec,Vec,Vec,void*);
/**
 * @brief The Jacobian evaluation routine for the dynamics simulation
 * @param [unknown] TS ts
 * @param [unknown] PetscReal t  
 * @param [unknown] Vec X 
 * @param [unknown] Vec Xdot 
 * @param [unknown] PetscReal a 
 * @param [unknown] Mat J
 * @param [unknown] Mat Jpre
 * @param [unknown] void* ctx
 * Notes: See PETSc routine TSSetIJacobian()
 */
extern PetscErrorCode DYNIJacobian(TS,PetscReal,Vec,Vec,PetscReal,Mat,Mat,void*);
/**
 * @brief Initializes the network and the dynamic state variables in the solution vector X
 * @param [in] DYN DYN - the dynamics simulation application object
 * @param [out] Vec X - the solution vector  
 */
extern PetscErrorCode DYNSetInitialConditions(DYN,Vec);
/**
 * @brief Computes the sensitivity of the dynamic constraints w.r.t. the parameters
 * @param [in] DYN DYN - the dynamics simulation application object  
 */
extern PetscErrorCode DYNComputeSensiP(DYN);
#ifdef FWDSA
/**
 * @brief Computes the Jacobian of the DYN equations w.r.t. parameters for the forward model
 * @param [in] TS ts - the TS solver 
 * @param [in] PetscReal t  - the current time
 * @param [in] Vec x  - the solution vector 
 * @param [out] Vec* jacP - Jacobian of the DAE equations w.r.t. parameters  
 * @param [in] void* ctx - application context (dyn)
 */
extern PetscErrorCode DYNComputeForwardJacobianP(TS,PetscReal,Vec,Vec*,void*);
#endif
/**
 * @brief Residual function for SNES to solve the algebraic part when an event is triggered
 * @param [unknown] SNES snes
 * @param [unknown] Vec X
 * @param [unknown] Vec F 
 * @param [unknown] void* ctx
 */
extern PetscErrorCode DYNEventPostSolveAlgebraicSensi(SNES,Vec,Vec,void*);
/**
 * @brief Jacobian routine for the post event algebraic solve
 * @param [unknown] SNES snes
 * @param [unknown] Vec X
 * @param [unknown] Mat J 
 * @param [unknown] Mat Jpre 
 * @param [unknown] void* ctx
 */
extern PetscErrorCode DYNEventPostSolveAlgebraicJacobian(SNES,Vec,Mat,Mat,void*);
#endif

extern PetscErrorCode DYNPostStage(TS,PetscReal,PetscInt,Vec*);
extern PetscErrorCode DYNPostEvaluate(TS);
