#include <exago_config.h>

#if defined(EXAGO_HAVE_RAJA)

#ifndef _PBPOLRAJAHIOP_H
#define _PBPOLRAJAHIOP_H

#include <opflow.h>

struct _p_FormPBPOLRAJAHIOP{};
typedef struct _p_FormPBPOLRAJAHIOP *PBPOLRAJAHIOP;

extern PetscErrorCode OPFLOWSetVariableBounds_PBPOLRAJAHIOP(OPFLOW,Vec,Vec);
extern PetscErrorCode OPFLOWSetVariableBoundsArray_PBPOLRAJAHIOP(OPFLOW,double*,double*);
extern PetscErrorCode OPFLOWSetConstraintBounds_PBPOLRAJAHIOP(OPFLOW,Vec,Vec);
extern PetscErrorCode OPFLOWSetConstraintBoundsArray_PBPOLRAJAHIOP(OPFLOW,double*,double*);
extern PetscErrorCode OPFLOWComputeEqualityConstraints_PBPOLRAJAHIOP(OPFLOW,Vec,Vec);
extern PetscErrorCode OPFLOWComputeInequalityConstraints_PBPOLRAJAHIOP(OPFLOW,Vec,Vec);
extern PetscErrorCode OPFLOWSetInitialGuessArray_PBPOLRAJAHIOP(OPFLOW,double*);
extern PetscErrorCode OPFLOWComputeEqualityConstraintsArray_PBPOLRAJAHIOP(OPFLOW,const double*,double*);
extern PetscErrorCode OPFLOWComputeInequalityConstraintsArray_PBPOLRAJAHIOP(OPFLOW,const double*,double*);
extern PetscErrorCode OPFLOWComputeObjective_PBPOLRAJAHIOP(OPFLOW,Vec,PetscScalar*);
extern PetscErrorCode OPFLOWComputeGradient_PBPOLRAJAHIOP(OPFLOW,Vec,Vec);
extern PetscErrorCode OPFLOWComputeObjectiveArray_PBPOLRAJAHIOP(OPFLOW,const double*,double*);
extern PetscErrorCode OPFLOWComputeGradientArray_PBPOLRAJAHIOP(OPFLOW,const double*,double*);
extern PetscErrorCode OPFLOWComputeEqualityConstraintJacobian_PBPOLRAJAHIOP(OPFLOW,Vec,Mat);
extern PetscErrorCode OPFLOWComputeInequalityConstraintJacobian_PBPOLRAJAHIOP(OPFLOW,Vec,Mat);
extern PetscErrorCode OPFLOWComputeObjectiveHessian_PBPOLRAJAHIOP(OPFLOW,Vec,Mat);
extern PetscErrorCode OPFLOWComputeEqualityConstraintsHessian_PBPOLRAJAHIOP(OPFLOW,Vec,Vec,Mat);
extern PetscErrorCode OPFLOWComputeInequalityConstraintsHessian_PBPOLRAJAHIOP(OPFLOW,Vec,Vec,Mat);
extern PetscErrorCode OPFLOWComputeSparseJacobian_PBPOLRAJAHIOP(OPFLOW,int*,int*,double*);
extern PetscErrorCode OPFLOWComputeSparseHessian_PBPOLRAJAHIOP(OPFLOW,const double*,int*,int*,double*);
extern PetscErrorCode OPFLOWComputeDenseEqualityConstraintJacobian_PBPOLRAJAHIOP(OPFLOW,const double*,double**);
extern PetscErrorCode OPFLOWComputeDenseInequalityConstraintJacobian_PBPOLRAJAHIOP(OPFLOW,const double*,double**);
extern PetscErrorCode OPFLOWComputeDenseHessian_PBPOLRAJAHIOP(OPFLOW,const double*,const double*,double**);

#endif
#endif
