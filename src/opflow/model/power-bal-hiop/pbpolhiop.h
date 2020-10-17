#include <exago_config.h>

#if defined(EXAGO_HAVE_HIOP)

#ifndef _PBPOLHIOP_H
#define _PBPOLHIOP_H

#include <opflow.h>
#include "../power-bal-polar2/pbpol2.h"

typedef struct _p_FormPBPOLHIOP *PBPOLHIOP;

struct _p_FormPBPOLHIOP{
  GENParams genparams;
  LOADParams loadparams;
  LINEParams lineparams;
  BUSParams  busparams;
};


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
extern PetscErrorCode OPFLOWComputeDenseEqualityConstraintJacobian_PBPOLHIOP(OPFLOW,const double*,double**);
extern PetscErrorCode OPFLOWComputeDenseInequalityConstraintJacobian_PBPOLHIOP(OPFLOW,const double*,double**);
extern PetscErrorCode OPFLOWComputeDenseHessian_PBPOLHIOP(OPFLOW,const double*,const double*,double**);

#endif
#endif
