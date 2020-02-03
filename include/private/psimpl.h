/**
 * @file psimpl.h
 * @brief Private header file containing the data for the power system stored into BUS, LINE, GEN, and LOAD data structures. It also contains its private API declarations
 */
#ifndef __PSIMPL_H
#define __PSIMPL_H

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
* @brief private bus data struct 
*/
struct _p_PSBUS{
  PetscInt      bus_i; /**< Integer bus number */
  char          i[20]; /**< Bus Number */
  char          name[20]; /**< Bus Name */
  PetscScalar   basekV; /**< Bus Base kV */
  PetscInt      ide; /**< Bus type code */
  PetscScalar   gl; /**< Active component of shunt admittance to ground */
  PetscScalar   bl; /**< Reactive component of shunt admittance to ground */
  PetscInt      area; /**< Area number */
  PetscInt      zone; /**< Zone number */
  PetscScalar   vm; /**< Bus voltage magnitude; in pu */
  PetscScalar   va; /**< Bus voltage phase angle */
  PetscInt      owner; /**< Owner number */
  PetscScalar   Vmax; /**< Max. voltage magnitude */
  PetscScalar   Vmin; /**< Min. voltage magnitude */
  PetscScalar   nvhi; /**< Normal voltage magnitude high limit */
  PetscScalar   nvlo; /**< Normal voltage magnitude low limit */
  PetscScalar   evhi; /**< Emergency voltage magnitude high limit */
  PetscScalar   evlo; /**< Emergency voltage magnitude low limit */
  PetscInt      internal_i; /**< Internal Bus Number */
  PetscInt      ngen; /**< Number of generators incident at this bus */
  PetscInt      ngenON; /**< Number of active generators at this bus */
  PetscInt      gidx[NGEN_AT_BUS_MAX]; /**< list of inndices for accessing the generator data in GEN structure */
  PetscInt      nload; /**< Number of loads incident at this bus */
  PetscInt      lidx[NLOAD_AT_BUS_MAX];
  PetscInt      nshunt; /**< Number of fixed shunts at this bus */
  PetscScalar   qrange; /**< Sum(Qmax_i - Qmin_i) for all ON generators at this bus */
  PetscScalar   qmintot; /**< Sum(Qmin_i) for all ON generators at this bus */
  PetscScalar   Pgtot; /**< Sum(Pg) for all ON generators at this bus, needed to split Pg at slack and Qg at all generator buses */
  PetscScalar   MVAbasetot; /**< Sum(MVAbase) for all ON generators at this bus, needed to split Pg at slack; Qg at all generator buses */

 
  PetscInt      nconnlines;  /**< Number of connected edges to this bus */
  PSLINE        connlines[MAXCONNLINES];   /**< Array of connected edges to this bus */

  
  PetscBool isghost;  /**< Is the bus ghosted? (Owned by other processor) */

  PSGEN     gens[NGEN_AT_BUS_MAX];  /**< Array of included generators to this bus */
  PSLOAD    loads[NLOAD_AT_BUS_MAX];  /**< Array of connected loads to this bus */

  
  PetscInt  startloc; /**< Starting location for the variables for this bus in the application residual "local" array.
     The local array also contains ghosted values. startloc is typically used to access values in
     the local vector, while startlocglob is used for setting the values in the matrix */
  
  PetscInt  startlocglob; /**< Starting location for the variables for this bus in the application "global" residual array */
};

/**
* @brief private load data struct
*/
struct _p_PSLOAD{
  PetscInt      bus_i; /**< Bus number */
  char          i[20]; /**< Bus Number or extended bus name*/
  char          id[3]; /**< Load identifier, in case of multiple loads. 1 by default */
  PetscInt      status; /**< Load status */
  PetscInt      area; /**< Area to which load is assigned */
  PetscInt      zone; /**< Zone to which load is assigned */
  PetscScalar   pl; /**< Active power component of constant MVA load */
  PetscScalar   ql; /**< Reactive power component of constant MVA load */
  PetscScalar   ip; /**< Active power component of constant current load: MW pu V */
  PetscScalar   iq; /**< Reactive power component of constant current load: Mvar pu V */
  PetscScalar   yp; /**< Active power component of constant admittance load: MW pu V */
  PetscScalar   yq; /**< Reactive power component of constant admittance load: Mvar pu V */
  PetscInt      owner; /**< Owner number */
  PetscInt      internal_i; /**< Internal Bus Number */
  PetscInt      scale;  /* Load scaling flag of one for a scalable load and zero for a fixed load */
  PetscInt      intrpt; /**< Interruptible load flag of one for an interruptible load for zero for a non interruptible load */
};

/**
* @brief private generator data struct 
* 20+ columns 
* 20, USING ONLY 1 OWNER's WORTH OF DATA. COME BACK TO THIS LATER, if necessary
*/
struct _p_PSGEN{
  PetscInt      bus_i;
  char          i[20]; /**< Bus Number or extended bus name*/
  char          id[3]; /**< Generator identifier, in case of multiple generators at same bus. 1 by default */
  PetscScalar   pg; /**< Generator active power output */
  PetscScalar   qg; /**< Generator reactive power output */
  PetscScalar   qt; /**< Maximum reactive power output: Mvar */
  PetscScalar   qb; /**< Minimum reactive power output: Mvar */
  PetscScalar   vs; /**< Regulated voltage setpoint: pu */
  PetscInt      ireg; /**< Remote bus number/identifier */
  PetscScalar   mbase; /**< MVA base of the machine */
  PetscScalar   zr; /**< Complex machine impedance ZSOURCE in pu on mbase */
  PetscScalar   zx; /**< ----------------------"------------------------- */
  PetscScalar   rt; /**< Step-up transformer impedance XTRAN in pu on mbase */
  PetscScalar   xt; /**< -----------------------"-------------------------- */
  PetscScalar   gtap; /**< Step-up transformer turns ratio */
  PetscInt      status; /**< Machine status */
  /**< Initial status of the generator. This is set once at the start of the simulation.
     During the simulation, the generator status can be updated
  */
  PetscInt      initial_status; 
  PetscScalar   rmpct; /**< Mvar % required to hold voltage at remote bus */
  PetscScalar   pt; /**< Gen max active power output: MW */
  PetscScalar   pb; /**< Gen min active power output: MW */
  PetscInt      o1; /**< Owner number */
  PetscScalar   f1; /**< Fraction of ownership */
  PetscInt      internal_i; /**< Internal Bus Number */
  PetscScalar   scale_gen;
  /* Generator cost data */
  PetscInt      cost_model; /**< generator cost model */
  PetscScalar   cost_startup; /**< generator startup cost */
  PetscScalar   cost_shutdown; /**< generator shutdown cost */

  /* Quadratic cost function is alpha*Pg^2 + beta*Pg + gamma */
  PetscInt      cost_ncoeffs; /**< Number of cost coefficients */
  PetscScalar   cost_gamma;
  PetscScalar   cost_beta;
  PetscScalar   cost_alpha;

  PetscBool     hasexc; /**< Flag for checking if this generator has an exciter */
  PetscBool     hasturbgov; /**< Flag for checking if this generator has a turbine-governor */
  PetscBool     hasstab;    /**< Flag for checking if this generator has a stabilizer */
};

/**
* @brief private line data structure
* 17+ columns 
*/
struct _p_PSLINE{
  PetscInt      fbus;
  PetscInt      tbus;
  char          i[20]; /**< Bus Number or extended bus name*/
  char          j[20]; /**< Bus Number or extended bus name*/
  char          ckt[20]; /**< Circuit identifier. 1 by default */
  PetscScalar   r; /**< Branch resistance: pu */
  PetscScalar   x; /**< Branch reactance: pu */
  PetscScalar   b; /**< Branch charging susceptance: pu */
  PetscScalar   rateA; /**< rate A in MVA */
  PetscScalar   rateB; /**< rate B in MVA */
  PetscScalar   rateC; /**< rate C in MVA */
  PetscScalar   tapratio;
  PetscScalar   phaseshift;
  PetscScalar   gi; /**< Complex admittance at 'i' end: pu */
  PetscScalar   bi; /**< Complex admittance at 'i' end: pu */
  PetscScalar   gj; /**< Complex admittance at 'j' end: pu */
  PetscScalar   bj; /**< Complex admittance at 'j' end: pu */
  PetscInt      status; /**< Service status */
  PetscInt      met;  /* Metered end flag; <= 1 to designate bus I as the metered end >= 2 to designate bus J as the metered en */
  PetscScalar   length; /**< Line length */
  PetscInt      o1; /**< Owner number */
  PetscScalar   f1; /**< Fraction of ownership */
  PetscScalar   sbase12;
  PetscInt      internal_i; /**< Internal From Bus Number */
  PetscInt      internal_j; /**< Internal To Bus Number */
  PetscScalar   yff[2],yft[2],ytf[2],ytt[2]; /**< [G,B] */
  PetscScalar   pf,qf,pt,qt; /**< Real and reactive power flows from and to ends */
  PetscBool     reversed_ends; /* Reversed line end (bus numbers are swapped) */

  PSBUS  connbuses[2]; /**< From and to buses */

  
  PetscInt startloc; /**< Starting location for the variables for this line in the local vector */

  
  PetscInt startlocglob; /**< Starting location for the variables for this line in the global vector */
};

typedef struct {
  PetscInt nv; /* Number of vertices in connected set */
  PetscInt *v; /* Vertices in the set */
  PetscInt blackout; /* No generation in this island, all buses isolated */
}PSConnComp;

typedef struct {
  PetscInt rank;
  PetscInt nc;
  PetscInt *c;
}PSConnCompgroup;

typedef struct {
    PetscInt        nc;
    PSConnCompgroup **cg; /* Pointers to connected groups */
}PSConngroupi;

typedef struct {
  PetscInt     n;   /* Number of connected component groups */
  PSConngroupi *ci;
}PSConngroup;
/**
* @brief private base power system data structure
*/
struct _p_PS {
  PetscScalar MVAbase;                 /* System base MVA */
  PetscInt    Nbus,Ngen,Nbranch,Nload; /* global # of buses,gens,branches, and loads (includes elements which are
                                          out of service */
  PetscInt    nbus,ngen,nbranch,nload; /* local # of buses,gens,branches,and loads */
  PetscInt    NgenON,ngenON;           /* Number of active generators */
  PetscInt    Nref,nref;               /* Number of reference buses */
  PSBUS       bus;
  PSLOAD      load;
  PSGEN       gen;
  PSLINE      line;
 
  COMM        comm;           /* communicator context */
  PetscInt    refct;          /* Reference count to keep track on how many objects are sharing PS */
  PetscInt    *busext2intmap; /* Maps external bus numbers to internal bus numbers */
  PetscInt    maxbusnum;      /* Max. bus number -- used for allocating busext2intmap */

  void*       app;            /* the application using this ps */
  PSApp       appname;        /* the application name using this ps */

  PetscInt    ndiff;          /* Number of differential equations.. only used for applications involving differential eqs. */
  
  PetscInt    compkey[10];    /* keys for components */
  DM          networkdm;      /* DM for managing the network */

  char        net_file_name[1024]; /* Network file name */
  PetscInt    nconncomp;           /* Number of connected components */
  PSConnComp  conncomp[MAXCONNECTEDCOMPS]; /* List of connected components */

  PetscBool setupcalled; /* Is setup called on PS? */
};

extern PetscErrorCode PSCheckTopology(PS);
extern PetscErrorCode PSGetLine(PS,PetscInt,PetscInt,const char*,PSLINE*);
extern PetscErrorCode PSGetGen(PS,PetscInt,const char*,PSGEN*);
extern PetscErrorCode PSSetGenStatus(PS,PetscInt,const char*,PetscInt);
extern PetscErrorCode PSSetLineStatus(PS,PetscInt,PetscInt,const char*,PetscInt);
extern PetscErrorCode PSIslandCheckandSetRefBus(PS,PetscInt);
extern PetscErrorCode PSConnCompDestroy(PS);
extern PetscErrorCode PSSetEdgeandBusStartLoc(PS);
#endif

