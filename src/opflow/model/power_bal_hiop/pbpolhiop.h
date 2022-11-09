
#include <exago_config.h>

#if defined(EXAGO_ENABLE_HIOP)

#ifndef _PBPOLHIOP_H
#define _PBPOLHIOP_H

#include <opflow.h>

typedef struct {
  int nbus;        /* Number of buses */
  int *isref;      /* isref[i] = 1 if bus is reference bus */
  int *isisolated; /* isisolated[i] = 1 if bus is isolated bus */
  int *ispvpq;     /* For all other buses */
  double *vmin;    /* min. voltage magnitude limit */
  double *vmax;    /* max. voltage magnitude limit */
  double *va;      /* bus angle (from file only used in bounds) */
  double *vm;      /* bus voltage magnitude (from file only used in bounds) */
  double *gl;      /* bus shunt (conductance) */
  double *bl;      /* bus shunt (suspectance) */
  double *powerimbalance_penalty; /* Penalty for power imbalance */
  int *xidx;     /* starting locations for bus variables in X vector */
  int *gidx;     /* starting locations for bus balance equations in constraint
                    vector */
  int *xidxpimb; /* starting locations for power imbalance bus variables */

  /* The following members are only used with HIOP */
  int *jacsp_idx;  /* Location number in the sparse Jacobian for Pimb */
  int *jacsq_idx;  /* Location number in the sparse Jacobian for Qimb */
  int *hesssp_idx; /* Location number in the Hessian */
} BUSParams;

typedef struct {
  int ngenON;         /* Number of generators with STATUS ON */
  double *cost_alpha; /* generator cost coefficients */
  double *cost_beta;  /* generator cost coefficients */
  double *cost_gamma; /* generator cost coefficients */
  double *pt;         /* min. active power gen. limits */
  double *pb;         /* max. active power gen. limits */
  double *qt;         /* min. reactive power gen. limits */
  double *qb;         /* max. reactive power gen. limits */
  double *pgs;        /* real power output setpoint */
  int *xidx;          /* starting locations in X vector */
  int *isrenewable;   /* Is a renewable generator? */
  int *
      gidxbus; /* starting locations in constraint vector for bus constraints */
  int *geqidxgen;    /* starting locations in equality constraint vector for gen
                        constraints */
  int *gineqidxgen;  /* starting locations in inequality constraint vector for
                        gen constraints */
  int *gbineqidxgen; /* Starting location to insert contribution to inequality
                        constraint bound */
  /* The following members are only used with HIOP */
  int *eqjacspbus_idx; /* Location number in the bus equality constraints sparse
                          Jacobian for Pg */
  int *eqjacsqbus_idx; /* Location number in the bus equality constraints sparse
                          Jacobian for Qg */
  int *eqjacspgen_idx; /* Location number in the gen equality constraints sparse
                          Jacobian for Pg */
  int *ineqjacspgen_idx; /* Location number in the bus equality constraints
                            sparse Jacobian for Pg */
  int *hesssp_idx;       /* Location number in the Hessian */
} GENParams;

typedef struct {
  int nload;                /* Number of loads */
  double *pl;               /* active power demand */
  double *ql;               /* reactive power demand */
  double *loadloss_penalty; /* Penalty for load loss */
  int *xidx;                /* starting location in X vector */
  int *gidx;                /* starting location in constraint vector */

  /* The following members are only used with HIOP */
  int *jacsp_idx;  /* Location number in the sparse Jacobian for delPload */
  int *jacsq_idx;  /* Location number in the sparse Jacobian for delQload */
  int *hesssp_idx; /* Location number in the Hessian */
} LOADParams;

typedef struct {
  int nlineON;   /* Number of active lines (STATUS = 1) */
  int nlinelim;  /* Active lines + limits */
  double *Gff;   /* From side self conductance */
  double *Bff;   /* From side self susceptance */
  double *Gft;   /* From-to transfer conductance */
  double *Bft;   /* From-to transfer susceptance */
  double *Gtf;   /* To-from transfer conductance */
  double *Btf;   /* To-from transfer susceptance */
  double *Gtt;   /* To side self conductance */
  double *Btt;   /* To side self susceptance */
  double *rateA; /* Line MVA rating A (normal operation) */
  int *xidxf;    /* Starting locatin of from bus voltage variables */
  int *xidxt;    /* Starting location of to bus voltage variables */
  int *geqidxf;  /* Starting location of from side to insert equality constraint
                    contribution in constraints vector */
  int *geqidxt;  /* Starting location of to side to insert equality constraint
                    contribution in constraints vector */
  int *gineqidx; /* Starting location to insert contribution to inequality
                    constraint */
  int *gbineqidx;  /* Starting location to insert contribution to inequality
                      constraint bound */
  int *linelimidx; /* Indices for subset of lines that have finite limits */
} LINEParams;

typedef struct _p_FormPBPOLHIOP *PBPOLHIOP;

struct _p_FormPBPOLHIOP {
  GENParams genparams;
  LOADParams loadparams;
  LINEParams lineparams;
  BUSParams busparams;
};

extern PetscErrorCode CreateBusParams(OPFLOW, BUSParams *);
extern PetscErrorCode CreateLineParams(OPFLOW, LINEParams *);
extern PetscErrorCode CreateLoadParams(OPFLOW, LOADParams *);
extern PetscErrorCode CreateGenParams(OPFLOW, GENParams *);

extern PetscErrorCode DestroyBusParams(OPFLOW, BUSParams *);
extern PetscErrorCode DestroyLineParams(OPFLOW, LINEParams *);
extern PetscErrorCode DestroyLoadParams(OPFLOW, LOADParams *);
extern PetscErrorCode DestroyGenParams(OPFLOW, GENParams *);

extern PetscErrorCode OPFLOWSetInitialGuessArray_PBPOLHIOP(OPFLOW, double *);
extern PetscErrorCode OPFLOWSetVariableBoundsArray_PBPOLHIOP(OPFLOW, double *,
                                                             double *);
extern PetscErrorCode OPFLOWSetConstraintBoundsArray_PBPOLHIOP(OPFLOW, double *,
                                                               double *);

extern PetscErrorCode
OPFLOWComputeEqualityConstraintsArray_PBPOLHIOP(OPFLOW, const double *,
                                                double *);
extern PetscErrorCode
OPFLOWComputeInequalityConstraintsArray_PBPOLHIOP(OPFLOW, const double *,
                                                  double *);

extern PetscErrorCode
OPFLOWComputeObjectiveArray_PBPOLHIOP(OPFLOW, const double *, double *);
extern PetscErrorCode
OPFLOWComputeGradientArray_PBPOLHIOP(OPFLOW, const double *, double *);
extern PetscErrorCode
OPFLOWComputeSparseEqualityConstraintJacobian_PBPOLHIOP(OPFLOW, const double *,
                                                        int *, int *, double *);
extern PetscErrorCode OPFLOWComputeSparseInequalityConstraintJacobian_PBPOLHIOP(
    OPFLOW, const double *, int *, int *, double *);
extern PetscErrorCode
OPFLOWComputeSparseHessian_PBPOLHIOP(OPFLOW, const double *, const double *,
                                     int *, int *, double *);
extern PetscErrorCode
OPFLOWComputeDenseEqualityConstraintJacobian_PBPOLHIOP(OPFLOW, const double *,
                                                       double *);
extern PetscErrorCode
OPFLOWComputeDenseInequalityConstraintJacobian_PBPOLHIOP(OPFLOW, const double *,
                                                         double *);
extern PetscErrorCode OPFLOWComputeDenseHessian_PBPOLHIOP(OPFLOW,
                                                          const double *,
                                                          const double *,
                                                          double *);

extern PetscErrorCode
OPFLOWSolutionCallback_PBPOLHIOP(OPFLOW, const double *, const double *,
                                 const double *, const double *, const double *,
                                 double obj_value);

#endif
#endif
