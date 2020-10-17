#ifndef _PBPOLRAJA_H
#define _PBPOLRAJA_H

#include <opflow.h>
#include <exago_config.h>

#if defined(EXAGO_HAVE_RAJA)

// Optimal power flow model class declaration
struct _p_FormPBPOLRAJA{};
typedef struct _p_FormPBPOLRAJA *PBPOLRAJA;

// Optimal power flow class public methods
extern PetscErrorCode OPFLOWSetVariableBounds_PBPOLRAJA(OPFLOW,Vec,Vec);
extern PetscErrorCode OPFLOWSetConstraintBounds_PBPOLRAJA(OPFLOW,Vec,Vec);
extern PetscErrorCode OPFLOWComputeEqualityConstraints_PBPOLRAJA(OPFLOW,Vec,Vec);
extern PetscErrorCode OPFLOWComputeInequalityConstraints_PBPOLRAJA(OPFLOW,Vec,Vec);
extern PetscErrorCode OPFLOWComputeObjective_PBPOLRAJA(OPFLOW,Vec,PetscScalar*);
extern PetscErrorCode OPFLOWComputeGradient_PBPOLRAJA(OPFLOW,Vec,Vec);
extern PetscErrorCode OPFLOWComputeEqualityConstraintJacobian_PBPOLRAJA(OPFLOW,Vec,Mat);
extern PetscErrorCode OPFLOWComputeInequalityConstraintJacobian_PBPOLRAJA(OPFLOW,Vec,Mat);
extern PetscErrorCode OPFLOWComputeObjectiveHessian_PBPOLRAJA(OPFLOW,Vec,Mat);
extern PetscErrorCode OPFLOWComputeEqualityConstraintsHessian_PBPOLRAJA(OPFLOW,Vec,Vec,Mat);
extern PetscErrorCode OPFLOWComputeInequalityConstraintsHessian_PBPOLRAJA(OPFLOW,Vec,Vec,Mat);

#endif

#endif
