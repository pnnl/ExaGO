
#include <exago_config.h>

#if defined(EXAGO_HAVE_HIOP)

#ifndef _PBPOLHIOP_H
#define _PBPOLHIOP_H

#include <opflow.h>

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

typedef struct _p_FormPBPOLHIOP *PBPOLHIOP;

struct _p_FormPBPOLHIOP{
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

extern PetscErrorCode OPFLOWSetVariableBounds_PBPOLHIOP(OPFLOW,Vec,Vec);
extern PetscErrorCode OPFLOWSetVariableBoundsArray_PBPOLHIOP(OPFLOW,double*,double*);
extern PetscErrorCode OPFLOWSetConstraintBounds_PBPOLHIOP(OPFLOW,Vec,Vec);
extern PetscErrorCode OPFLOWSetConstraintBoundsArray_PBPOLHIOP(OPFLOW,double*,double*);
extern PetscErrorCode OPFLOWComputeEqualityConstraints_PBPOLHIOP(OPFLOW,Vec,Vec);
extern PetscErrorCode OPFLOWComputeInequalityConstraints_PBPOLHIOP(OPFLOW,Vec,Vec);
extern PetscErrorCode OPFLOWComputeEqualityConstraintsArray_PBPOLHIOP(OPFLOW,const double*,double*);
extern PetscErrorCode OPFLOWComputeInequalityConstraintsArray_PBPOLHIOP(OPFLOW,const double*,double*);
extern PetscErrorCode OPFLOWComputeObjective_PBPOLHIOP(OPFLOW,Vec,PetscScalar*);
extern PetscErrorCode OPFLOWComputeGradient_PBPOLHIOP(OPFLOW,Vec,Vec);
extern PetscErrorCode OPFLOWComputeObjectiveArray_PBPOLHIOP(OPFLOW,const double*,double*);
extern PetscErrorCode OPFLOWComputeGradientArray_PBPOLHIOP(OPFLOW,const double*,double*);
extern PetscErrorCode OPFLOWComputeEqualityConstraintJacobian_PBPOLHIOP(OPFLOW,Vec,Mat);
extern PetscErrorCode OPFLOWComputeInequalityConstraintJacobian_PBPOLHIOP(OPFLOW,Vec,Mat);
extern PetscErrorCode OPFLOWComputeObjectiveHessian_PBPOLHIOP(OPFLOW,Vec,Mat);
extern PetscErrorCode OPFLOWComputeEqualityConstraintsHessian_PBPOLHIOP(OPFLOW,Vec,Vec,Mat);
extern PetscErrorCode OPFLOWComputeInequalityConstraintsHessian_PBPOLHIOP(OPFLOW,Vec,Vec,Mat);
extern PetscErrorCode OPFLOWComputeSparseJacobian_PBPOLHIOP(OPFLOW,int*,int*,double*);
extern PetscErrorCode OPFLOWComputeSparseHessian_PBPOLHIOP(OPFLOW,const double*,int*,int*,double*);
extern PetscErrorCode OPFLOWComputeDenseEqualityConstraintJacobian_PBPOLHIOP(OPFLOW,const double*,double*);
extern PetscErrorCode OPFLOWComputeDenseInequalityConstraintJacobian_PBPOLHIOP(OPFLOW,const double*,double*);
extern PetscErrorCode OPFLOWComputeDenseHessian_PBPOLHIOP(OPFLOW,const double*,const double*,double*);

#endif
#endif
