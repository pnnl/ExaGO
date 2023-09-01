/**
 * @file psimpl.h
 * @brief Private header file containing the data for the power system stored
 * into BUS, LINE, GEN, and LOAD data structures. It also contains its private
 * API declarations
 */
#ifndef __PSIMPL_H
#define __PSIMPL_H

#include <private/contingencylist.h>
#include <private/scenariolist.h>
#include <ps.h>

/**
 * @brief Maximum number of lines allowed to be connected to a bus
 */
#define MAXCONNLINES 64

/**
 * @brief Maximum number of connected components in the graph
 */
#define MAXCONNECTEDCOMPS 64

/**
 * @brief A bogus value for load loss cost
 */
#define BOGUSLOSSCOST -1234.0

/**
 * @brief private bus data struct
 */
struct _p_PSBUS {
  PetscInt bus_i;      /**< Integer bus number */
  char i[20];          /**< Bus Number */
  char name[20];       /**< Bus Name */
  PetscScalar basekV;  /**< Bus Base kV */
  PetscInt ide;        /**< Bus type code */
  PetscScalar gl;      /**< Active component of shunt admittance to ground */
  PetscScalar bl;      /**< Reactive component of shunt admittance to ground */
  PetscInt area;       /**< Area number */
  PetscInt zone;       /**< Zone number */
  PetscScalar vm;      /**< Bus voltage magnitude; in pu */
  PetscScalar va;      /**< Bus voltage phase angle */
  PetscInt owner;      /**< Owner number */
  PetscScalar Vmax;    /**< Max. voltage magnitude */
  PetscScalar Vmin;    /**< Min. voltage magnitude */
  PetscScalar nvhi;    /**< Normal voltage magnitude high limit */
  PetscScalar nvlo;    /**< Normal voltage magnitude low limit */
  PetscScalar evhi;    /**< Emergency voltage magnitude high limit */
  PetscScalar evlo;    /**< Emergency voltage magnitude low limit */
  PetscInt internal_i; /**< Internal Bus Number */
  PetscInt ngen;       /**< Number of generators incident at this bus */
  PetscInt ngenON;     /**< Number of active generators at this bus */
  PetscInt gidx[NGEN_AT_BUS_MAX]; /**< list of inndices for accessing the
                                     generator data in GEN structure */
  PetscInt nload;                 /**< Number of loads incident at this bus */
  PetscInt lidx[NLOAD_AT_BUS_MAX];
  PetscInt nshunt; /**< Number of fixed shunts at this bus */
  PetscScalar
      qrange; /**< Sum(Qmax_i - Qmin_i) for all ON generators at this bus */
  PetscScalar qmintot; /**< Sum(Qmin_i) for all ON generators at this bus */
  PetscScalar Pgtot;   /**< Sum(Pg) for all ON generators at this bus, needed to
                          split Pg at slack and Qg at all generator buses */
  PetscScalar
      MVAbasetot; /**< Sum(MVAbase) for all ON generators at this bus, needed to
                     split Pg at slack; Qg at all generator buses */

  PetscInt nconnlines;            /**< Number of connected edges to this bus */
  PSLINE connlines[MAXCONNLINES]; /**< Array of connected edges to this bus */

  PetscBool isghost; /**< Is the bus ghosted? (Owned by other processor) */

  PSGEN gens[NGEN_AT_BUS_MAX]; /**< Array of included generators to this bus */
  PSLOAD loads[NLOAD_AT_BUS_MAX]; /**< Array of connected loads to this bus */

  PetscScalar pimb; /* Real power imbalance (set from OPFLOW solution) */
  PetscScalar qimb; /* Reactive power imbalance (set from OPFLOW solution) */

  PetscScalar mult_pmis; /* Lagrange multiplier for Re(power/current balance)
                            mismatch */
  PetscScalar mult_qmis; /* Lagrange multiplier for Im(power/current balance)
                            mismatch */

  /* Variable, constraints sizes and locations */
  PetscInt startloc; /**< Starting location for the variables for this bus in
    the application residual "local" array. The local array also contains
    ghosted values. startloc is typically used to access values in the local
    vector, while startlocglob is used for setting the values in the matrix */
  PetscInt startlocglob; /**< Starting location for the variables for this bus
                            in the application "global" residual array */

  PetscInt nxV;     /* Number of variables for bus voltage */
  PetscInt nxshunt; /* Number of variables for bus shunt */
  PetscInt nxpimb;  /* Number of variables for power imbalance */

  PetscInt startxVloc;     /* Location for voltage variables */
  PetscInt startxshuntloc; /* Location for bus shunt variables */
  PetscInt startxpimbloc;  /* Location for power imbalance variables */

  PetscInt startxVlocglob;     /* Location for voltage variables */
  PetscInt startxshuntlocglob; /* Location for bus shunt variables */
  PetscInt startxpimblocglob;  /* Location for power imbalance variables */

  PetscInt nx; /* Total number of variables for the bus */

  PetscInt nconeq;   /* Number of equality constraints */
  PetscInt nconineq; /* Number of inequality constraints */

  PetscInt
      starteqloc; /* Starting location for equality constraints for this bus */
  PetscInt startineqloc; /* Starting location for inequality constraints for
                            this bus */
};

/**
 * @brief private load data struct
 */
struct _p_PSLOAD {
  PetscInt bus_i; /**< Bus number */
  char i[20];     /**< Bus Number or extended bus name*/
  char id[3]; /**< Load identifier, in case of multiple loads. 1 by default */
  PetscInt status; /**< Load status */
  PetscInt area;   /**< Area to which load is assigned */
  PetscInt zone;   /**< Zone to which load is assigned */
  PetscScalar pl;  /**< Active power component of constant MVA load */
  PetscScalar ql;  /**< Reactive power component of constant MVA load */
  PetscScalar
      ip; /**< Active power component of constant current load: MW pu V */
  PetscScalar
      iq; /**< Reactive power component of constant current load: Mvar pu V */
  PetscScalar
      yp; /**< Active power component of constant admittance load: MW pu V */
  PetscScalar yq; /**< Reactive power component of constant admittance load:
                     Mvar pu V */
  PetscInt owner; /**< Owner number */
  PetscInt internal_i; /**< Internal Bus Number */
  PetscInt scale; /* Load scaling flag of one for a scalable load and zero for a
                     fixed load */
  PetscInt intrpt; /**< Interruptible load flag of one for an interruptible load
                      for zero for a non interruptible load */

  PetscScalar pl_loss; /* Real power load loss */
  PetscScalar ql_loss; /* Reactive power load loss */

  PetscScalar loss_cost; /* Cost for the load..used in load shedding */
  PetscScalar loss_frac; /* Fraction of allowed load shed */

  /* Variable, constraint sizes and locations */

  PetscInt nxloadloss; /* Number of variables for load loss */

  PetscInt startxloadlossloc;     /* Location of load loss variables */
  PetscInt startxloadlosslocglob; /* Global Location of load loss variables */

  PetscInt nx; /* Total number of variables for the load */

  PetscInt nconeq;   /* Number of equality constraints */
  PetscInt nconineq; /* Number of inequality constraints */

  PetscInt startloc;     /* Starting location for variables */
  PetscInt startlocglob; /* Starting global location for variables */

  PetscInt starteqloc;   /* Starting location for equality constraints */
  PetscInt startineqloc; /* Starting location for inequality constraints */
};

/**
 * @brief private generator data struct
 * 20+ columns
 * 20, USING ONLY 1 OWNER's WORTH OF DATA. COME BACK TO THIS LATER, if necessary
 */
struct _p_PSGEN {
  PetscInt bus_i;
  char i[20]; /**< Bus Number or extended bus name*/
  char id[3]; /**< Generator identifier, in case of multiple generators at same
                 bus. 1 by default */
  PetscScalar pg;    /**< Generator active power output */
  PetscScalar qg;    /**< Generator reactive power output */
  PetscScalar qt;    /**< Maximum reactive power output: Mvar */
  PetscScalar qb;    /**< Minimum reactive power output: Mvar */
  PetscScalar vs;    /**< Regulated voltage setpoint: pu */
  PetscInt ireg;     /**< Remote bus number/identifier */
  PetscScalar mbase; /**< MVA base of the machine */
  PetscScalar zr;    /**< Complex machine impedance ZSOURCE in pu on mbase */
  PetscScalar zx;    /**< ----------------------"------------------------- */
  PetscScalar rt;    /**< Step-up transformer impedance XTRAN in pu on mbase */
  PetscScalar xt;    /**< -----------------------"-------------------------- */
  PetscScalar gtap;  /**< Step-up transformer turns ratio */
  PetscInt status;   /**< Machine status */
  /**< Initial status of the generator. This is set once at the start of the
     simulation. During the simulation, the generator status can be updated
  */
  PetscInt initial_status;
  PetscScalar rmpct;   /**< Mvar % required to hold voltage at remote bus */
  PetscScalar pt;      /**< Gen max active power output: MW */
  PetscScalar pb;      /**< Gen min active power output: MW */
  PetscInt o1;         /**< Owner number */
  PetscScalar f1;      /**< Fraction of ownership */
  PetscInt internal_i; /**< Internal Bus Number */
  PetscScalar scale_gen;
  /* Generator cost data */
  PetscInt cost_model;       /**< generator cost model */
  PetscScalar cost_startup;  /**< generator startup cost */
  PetscScalar cost_shutdown; /**< generator shutdown cost */

  /* Quadratic cost function is alpha*Pg^2 + beta*Pg + gamma */
  PetscInt cost_ncoeffs; /**< Number of cost coefficients */
  PetscScalar cost_gamma;
  PetscScalar cost_beta;
  PetscScalar cost_alpha;

  PetscScalar pc1;    /* lower real power output of PQ capability curve (pu) */
  PetscScalar pc2;    /* upper real power output of PQ capability curve (pu) */
  PetscScalar qc1min; /* minimum reactive power output at Pc1 (pu) */
  PetscScalar qc1max; /* maximum reactive power output at Pc1 (pu) */
  PetscScalar qc2min; /* minimum reactive power output at Pc2 (pu)  */
  PetscScalar qc2max; /* maximum reactive power output at Pc2 (pu) */
  PetscScalar
      ramp_rate_min; /* ramp rate or load following/AGC per minute (pu) */
  PetscScalar ramp_rate_10min;    /* ramp rate for 10 minute reserves (pu) */
  PetscScalar ramp_rate_30min;    /* ramp rate for 30 minute reserves (pu) */
  PetscScalar ramp_rate_min_mvar; /* ramp rate for reactive power (2 sec
                                     timescale) (pu) */
  PetscScalar apf;                /* area participation factor */

  PetscScalar pgs; /* Generator set-point power (used with scopflow, agc) */

  /* genfuel type - See GENFUEL_TYPES defined in constants.h */
  PetscInt genfuel_type;

  /* Variable, constraint sizes and locations */
  PetscInt nxpow;  /* Number of generator power variables (Pg, Qg)*/
  PetscInt nxpset; /* Number of generator set point variables (Pgset) */
  PetscInt nxpdev; /* Number of generator real power deviation (ramp up/down)
                      variables */

  PetscInt startxpowloc; /* Starting location for power variables */
  PetscInt
      startxpsetloc; /* Starting location of generator set-point variable */
  PetscInt startxpdevloc; /* Starting location for power deviation variables */

  PetscInt startxpowlocglob; /* Starting global location for power variables */
  PetscInt
      startxpsetlocglob; /* Starting location of generator set-point variable */
  PetscInt startxpdevlocglob; /* Starting global location for power deviation
                                 variables */

  PetscInt nx; /* Total number of variables for the gen */

  PetscInt nconeq;   /* Number of equality constraints */
  PetscInt nconineq; /* Number of inequality constraints */

  PetscInt startloc;     /* Starting location for variables */
  PetscInt starteqloc;   /* Starting location for equality constraints */
  PetscInt startineqloc; /* Starting location for inequality constraints */

  PetscBool isrenewable; /* Is this a renewable generator ? */
};

/**
 * @brief private line data structure
 * 17+ columns
 */
struct _p_PSLINE {
  PetscInt fbus;
  PetscInt tbus;
  char i[20];        /**< Bus Number or extended bus name*/
  char j[20];        /**< Bus Number or extended bus name*/
  char ckt[20];      /**< Circuit identifier. 1 by default */
  PetscScalar r;     /**< Branch resistance: pu */
  PetscScalar x;     /**< Branch reactance: pu */
  PetscScalar b;     /**< Branch charging susceptance: pu */
  PetscScalar rateA; /**< rate A in MVA */
  PetscScalar rateB; /**< rate B in MVA */
  PetscScalar rateC; /**< rate C in MVA */
  PetscScalar tapratio;
  PetscScalar phaseshift;
  PetscScalar gi;  /**< Complex admittance at 'i' end: pu */
  PetscScalar bi;  /**< Complex admittance at 'i' end: pu */
  PetscScalar gj;  /**< Complex admittance at 'j' end: pu */
  PetscScalar bj;  /**< Complex admittance at 'j' end: pu */
  PetscInt status; /**< Service status */
  PetscInt met; /* Metered end flag; <= 1 to designate bus I as the metered end
                   >= 2 to designate bus J as the metered en */
  PetscScalar length; /**< Line length */
  PetscInt o1;        /**< Owner number */
  PetscScalar f1;     /**< Fraction of ownership */
  PetscScalar sbase12;
  PetscInt internal_i;                        /**< Internal From Bus Number */
  PetscInt internal_j;                        /**< Internal To Bus Number */
  PetscScalar yff[2], yft[2], ytf[2], ytt[2]; /**< [G,B] */
  PetscScalar bdc,
      pshift; /**< line admittance and power injection from phase shiters */
  PetscScalar pf, qf, pt,
      qt;             /**< Real and reactive power flows from and to ends */
  PetscScalar sf, st; /**<Apparent power flows at the two ends */
  PetscBool reversed_ends; /* Reversed line end (bus numbers are swapped) */
  PetscScalar kvlevel; /* Kv level for lines, for transformers uses the HV side
                          voltage */

  PSBUS connbuses[2]; /**< From and to buses */

  /****** For DC lines **********/
  PetscBool isdcline; /**< Is line a DC line? */
  PetscScalar pmin; /**< lower limit on PF (MW flow at "from" end) */
  PetscScalar pmax; /**< upper limit on PF (MW flow at "from" end) */
  PetscScalar Vf;   /**<Voltage set-point at "from" bus (p.u.) */
  PetscScalar Vt;   /**<Voltage set-point at "to" bus (p.u.) */
  PetscScalar qminf; /**< lower limit on QF (MVAr flow at "from" end) */
  PetscScalar qmaxf; /**< upper limit on QF (MVAr flow at "from" end) */
  PetscScalar qmint; /**< lower limit on QT (MVAr flow at "to" end) */
  PetscScalar qmaxt; /**< upper limit on QT (MVAr flow at "to" end) */
  PetscScalar loss0,loss1; /* constant and linear term for loss function (loss = loss0 + loss1*PF) */
                           /* loss0 given in terms of MW, loss1 is dimensionless */
  PetscScalar
  mult_pmin,mult_pmax;
  PetscScalar
  mult_qminf,mult_qmaxf,mult_qmint,mult_qmaxt;
  /*******************/
  
  PetscScalar
      mult_sf; /* Lagrange multiplier for from bus injection (set by OPFLOW) */
  PetscScalar
      mult_st; /* Lagrange multiplier for to bus injection (set by OPFLOW) */

  // From and to end substations (used for visualization only)
  PSSUBST subst_from;
  PSSUBST subst_to;

  /* Variable and constraint sizes and locations */
  PetscInt startloc; /**< Starting location for the variables for this line in
                        the local vector */
  PetscInt startlocglob; /**< Starting location for the variables for this line
                            in the global vector */

  PetscInt nx; /* Number of variables for the line */

  PetscInt nconeq;   /* Number of equality constraints */
  PetscInt nconineq; /* Number of inequality constraints */

  PetscInt starteqloc;   /* Starting location for equality constraints */
  PetscInt startineqloc; /* Starting location for inequality constraints */
};

extern const char *const PSGENFuelTypes[];

typedef struct {
  PetscInt nv;       /* Number of vertices in connected set */
  PetscInt *v;       /* Vertices in the set */
  PetscInt blackout; /* No generation in this island, all buses isolated */
} PSConnComp;

typedef struct {
  PetscInt rank;
  PetscInt nc;
  PetscInt *c;
} PSConnCompgroup;

typedef struct {
  PetscInt nc;
  PSConnCompgroup **cg; /* Pointers to connected groups */
} PSConngroupi;

typedef struct {
  PetscInt n; /* Number of connected component groups */
  PSConngroupi *ci;
} PSConngroup;

struct _p_PSSUBST {
  PetscInt num;             /* Substation number */
  PetscInt intnum;          /* Internal number */
  char name[64];            /* Substation name */
  PetscScalar longlat[2];   /* longitude and latitude */
  PetscInt nbus;            /* Number of buses */
  PSBUS bus[20];            /* Pointers for buses */
  PetscInt nkvlevels;       /* Number of KV levels at this substation */
  PetscScalar kvlevels[10]; /* Substation KV levels */
};

/* Struct to save system summary stats */
typedef struct {
  PetscInt Nbus;             /* Number of buses */
  PetscInt Ngen;             /* Number of generators */
  PetscInt NgenON;           /* Number of committed generators */
  PetscInt Nline;            /* Number of lines (includes DC lines) */
  PetscInt NlineON;          /* Number of lines ON (includes DC lines) */
  PetscInt Ndcline;          /* Number of dc lines */
  PetscInt NdclineON;        /* Number of dc lines ON */

  PetscInt Nload;            /* Number of loads */
  PetscScalar total_pgencap; /* Total active generation capacity */

  PetscScalar total_genON[2];  /* Total generation online [MW, MVAr] */
  PetscScalar total_pgencapON; /* Total generation capacity online */

  PetscScalar total_load[2];     /* Total demand [MW, MVAr] */
  PetscScalar total_loadshed[2]; /* Total demand shed [MW, MVar] */
} PSSystemSummary;

/**
 * @brief private base power system data structure
 */
struct _p_PS {
  PetscScalar MVAbase; /* System base MVA */
  PetscInt Nbus, Ngen, Nline,
    Nload, Ndcline; /* global # of buses,gens,branches,loads, and dclines (includes elements
                which are out of service */
  PetscInt nbus, ngen, nline,
    nload,ndcline;                 /* local # of buses,gens,branches,and loads */
  PetscInt NlineON, nlineON; /* Number of active lines (includes DC lines) */
  PetscInt NdclineON, ndclineON; /* Number of active dc lines */
  PetscInt NgenON, ngenON;   /* Number of active generators */
  PetscInt Nref, nref;       /* Number of reference buses */
  /* Number of generator types */
  PetscInt ngencoal, ngenwind, ngensolar, ngenng, ngennuclear, ngenhydro,
      ngenundefined;
  /* Number of renewable generators (solar, wind) */
  PetscInt ngenrenew;
  /* Number of isolated buses */
  PetscInt nisolated_buses,Nisolated_buses;
  PetscInt *isolated_buses; /* Array to hold isolated buses */

  PSBUS bus;
  PSLOAD load;
  PSGEN gen;
  PSLINE line;

  COMM comm;      /* communicator context */
  PetscInt refct; /* Reference count to keep track on how many objects are
                     sharing PS */
  PetscInt
      *busext2intmap; /* Maps external bus numbers to internal bus numbers */
  PetscInt maxbusnum; /* Max. bus number -- used for allocating busext2intmap */

  void *app;     /* the application using this ps */
  PSApp appname; /* the application name using this ps */

  PetscInt ndiff; /* Number of differential equations.. only used for
                     applications involving differential eqs. */

  PetscInt compkey[10]; /* keys for components */
  DM networkdm;         /* DM for managing the network */

  char net_file_name[1024]; /* Network file name */

  char gic_file_name[1024];               /* GIC data file */
  PetscBool gic_file_set;                 /* Is GIC file set? */
  PetscInt nconncomp;                     /* Number of connected components */
  PSConnComp conncomp[MAXCONNECTEDCOMPS]; /* List of connected components */

  PetscBool opflow_converged; /* OPFLOW converged? */
  PetscReal opflowobj; /* OPFLOW objective value used when printing to file */

  /* Used by OPFLOW */
  /* System-level variables,constraints, and their location */
  PetscInt nx; /* Total number of variables for the bus */

  PetscInt nconeq;   /* Number of equality constraints */
  PetscInt nconineq; /* Number of inequality constraints */

  PetscInt startxloc;     /* Starting location for variables */
  PetscInt startxlocglob; /* Starting global location for variables */

  PetscInt
      starteqloc; /* Starting location for equality constraints for this bus */
  PetscInt startineqloc; /* Starting location for inequality constraints for
                            this bus */
  PetscInt nkvlevels;    /* Number of different kV levels */
  PetscScalar *kvlevels; /* kV levels */

  PSSUBST substations;
  PetscInt nsubstations;

  PSSystemSummary sys_info;

  PetscBool setupcalled; /* Is setup called on PS? */

  PetscLogDouble solve_real_time;
  PetscLogDouble solve_cpu_time;
};

extern PetscErrorCode PSCheckTopology(PS);
extern PetscErrorCode PSGetLine(PS, PetscInt, PetscInt, const char *, PSLINE *);
extern PetscErrorCode PSGetGen(PS, PetscInt, const char *, PSGEN *);
extern PetscErrorCode PSGetLoad(PS, PetscInt, const char *, PSLOAD *);
extern PetscErrorCode PSIslandCheckandSetRefBus(PS, PetscInt);
extern PetscErrorCode PSConnCompDestroy(PS);
extern PetscErrorCode PSSetEdgeandBusStartLoc(PS);
extern PetscErrorCode PSPrintSystemSummary(PS);
extern PetscErrorCode PSSaveSolution(PS, OutputFormat, const char[]);
extern PetscErrorCode PSApplyContingency(PS, struct _p_Contingency);
extern PetscErrorCode PSIncreaseReferenceCount(PS);
extern PetscErrorCode PSDescreaseReferenceCount(PS);
extern PetscErrorCode PSSetGenDispatchandStatus(PS, PetscInt, PetscInt,
                                                PetscInt, PetscScalar,
                                                PetscScalar);
extern PetscErrorCode PSApplyScenario(PS, Scenario);
extern PetscErrorCode PSGetKVLevels(PS, PetscInt *, const PetscScalar **);
extern PetscErrorCode PSReadGICData(PS);
#endif
