#include <exago_config.h>

#if defined(EXAGO_ENABLE_RAJA)

#ifndef _PBPOLRAJAHIOP_H
#define _PBPOLRAJAHIOP_H

#include <opflow.h>

struct _p_FormPBPOLRAJAHIOP{};
typedef struct _p_FormPBPOLRAJAHIOP *PBPOLRAJAHIOP;

extern PetscErrorCode OPFLOWSetVariableBoundsArray_PBPOLRAJAHIOP(OPFLOW,double*,double*);
extern PetscErrorCode OPFLOWSetConstraintBoundsArray_PBPOLRAJAHIOP(OPFLOW,double*,double*);
extern PetscErrorCode OPFLOWSetInitialGuessArray_PBPOLRAJAHIOP(OPFLOW,double*);
extern PetscErrorCode OPFLOWComputeEqualityConstraintsArray_PBPOLRAJAHIOP(OPFLOW,const double*,double*);
extern PetscErrorCode OPFLOWComputeInequalityConstraintsArray_PBPOLRAJAHIOP(OPFLOW,const double*,double*);
extern PetscErrorCode OPFLOWComputeObjectiveArray_PBPOLRAJAHIOP(OPFLOW,const double*,double*);
extern PetscErrorCode OPFLOWComputeGradientArray_PBPOLRAJAHIOP(OPFLOW,const double*,double*);
extern PetscErrorCode OPFLOWComputeSparseJacobian_PBPOLRAJAHIOP(OPFLOW,int*,int*,double*);
extern PetscErrorCode OPFLOWComputeSparseHessian_PBPOLRAJAHIOP(OPFLOW,const double*,int*,int*,double*);
extern PetscErrorCode OPFLOWComputeDenseEqualityConstraintJacobian_PBPOLRAJAHIOP(OPFLOW,const double*,double*);
extern PetscErrorCode OPFLOWComputeDenseInequalityConstraintJacobian_PBPOLRAJAHIOP(OPFLOW,const double*,double*);
extern PetscErrorCode OPFLOWComputeDenseHessian_PBPOLRAJAHIOP(OPFLOW,const double*,const double*,double*);

#endif
#endif
