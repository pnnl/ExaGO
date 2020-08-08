#include <opflow.h>

#ifndef _PBPOL2_H
#define _PBPOL2_H

typedef struct {
  int    nbus;    /* Number of buses */
  int    *isref;  /* isref[i] = 1 if bus is reference bus */
  int    *isisolated; /* isisolated[i] = 1 if bus is isolated bus */
  int    *ispvpq; /* For all other buses */
  double *vmin; /* min. voltage magnitude limit */
  double *vmax; /* max. voltage magnitude limit */
  double *va;   /* bus angle (from file only used in bounds) */
  double *vm;   /* bus voltage magnitude (from file only used in bounds) */
  double *gl;  /* bus shunt (conductance) */
  double *bl;  /* bus shunt (suspectance) */
  int    *xidx; /* starting locations for bus variables in X vector */
  int    *gidx;  /* starting locations for bus balance equations in constraint vector */
}BUSParams;

typedef struct {
  int    ngenON;       /* Number of generators with STATUS ON */
  double *cost_alpha;  /* generator cost coefficients */
  double *cost_beta;   /* generator cost coefficients */
  double *cost_gamma;  /* generator cost coefficients */
  double *pt;          /* min. active power gen. limits */
  double *pb;          /* max. active power gen. limits */
  double *qt;          /* min. reactive power gen. limits */
  double *qb;          /* max. reactive power gen. limits */
  int    *xidx;        /* starting locations in X vector */
  int    *gidx;         /* starting locations in constraint vector */

  /* The following members are only used with HIOP */
  int   *jacsp_idx; /* Location number in the sparse Jacobian for Pg */
  int   *jacsq_idx; /* Location number in the sparse Jacobian for Qg */
}GENParams;

typedef struct{
  int    nload; /* Number of loads */
  double *pl;   /* active power demand */
  double *ql;   /* reactive power demand */
  int    *xidx; /* starting location in X vector */
  int    *gidx;  /* starting location in constraint vector */
}LOADParams;

typedef struct{
  int     nlineON; /* Number of active lines (STATUS = 1) */
  int     nlinelim; /* Active lines + limits */
  double *Gff;  /* From side self conductance */
  double *Bff;  /* From side self susceptance */
  double *Gft;  /* From-to transfer conductance */
  double *Bft;  /* From-to transfer susceptance */
  double *Gtf;  /* To-from transfer conductance */
  double *Btf;  /* To-from transfer susceptance */
  double *Gtt;  /* To side self conductance */
  double *Btt;  /* To side self susceptance */
  double *rateA; /* Line MVA rating A (normal operation) */
  int    *xidxf; /* Starting locatin of from bus voltage variables */
  int    *xidxt; /* Starting location of to bus voltage variables */
  int    *geqidxf; /* Starting location of from side to insert equality constraint contribution in constraints vector */
  int    *geqidxt; /* Starting location of to side to insert equality constraint contribution in constraints vector */
  int    *gineqidx; /* Starting location to insert contribution to inequality constraint */
  int    *gbineqidx; /* Starting location to insert contribution to inequality constraint bound */
  int    *linelimidx; /* Indices for subset of lines that have finite limits */
}LINEParams;

typedef struct _p_FormPBPOL2 *PBPOL2;

struct _p_FormPBPOL2{
  GENParams genparams;
  LOADParams loadparams;
  LINEParams lineparams;
  BUSParams  busparams;
};

extern PetscErrorCode CreateBusParams(OPFLOW,BUSParams*);
extern PetscErrorCode CreateLineParams(OPFLOW,LINEParams*);
extern PetscErrorCode CreateLoadParams(OPFLOW,LOADParams*);
extern PetscErrorCode CreateGenParams(OPFLOW,GENParams*);

extern PetscErrorCode DestroyBusParams(BUSParams*);
extern PetscErrorCode DestroyLineParams(LINEParams*);
extern PetscErrorCode DestroyLoadParams(LOADParams*);
extern PetscErrorCode DestroyGenParams(GENParams*);

extern PetscErrorCode OPFLOWSetVariableBounds_PBPOL2(OPFLOW,Vec,Vec);
extern PetscErrorCode OPFLOWSetConstraintBounds_PBPOL2(OPFLOW,Vec,Vec);
extern PetscErrorCode OPFLOWComputeEqualityConstraints_PBPOL2(OPFLOW,Vec,Vec);
extern PetscErrorCode OPFLOWComputeInequalityConstraints_PBPOL2(OPFLOW,Vec,Vec);
extern PetscErrorCode OPFLOWComputeObjective_PBPOL2(OPFLOW,Vec,PetscScalar*);
extern PetscErrorCode OPFLOWComputeGradient_PBPOL2(OPFLOW,Vec,Vec);
extern PetscErrorCode OPFLOWComputeEqualityConstraintJacobian_PBPOL2(OPFLOW,Vec,Mat);
extern PetscErrorCode OPFLOWComputeInequalityConstraintJacobian_PBPOL2(OPFLOW,Vec,Mat);
extern PetscErrorCode OPFLOWComputeObjectiveHessian_PBPOL2(OPFLOW,Vec,Mat);
extern PetscErrorCode OPFLOWComputeEqualityConstraintsHessian_PBPOL2(OPFLOW,Vec,Vec,Mat);
extern PetscErrorCode OPFLOWComputeInequalityConstraintsHessian_PBPOL2(OPFLOW,Vec,Vec,Mat);

#endif
