#include <exago_config.h>

#if defined(EXAGO_ENABLE_RAJA)
#if defined(EXAGO_ENABLE_HIOP_SPARSE)

#ifndef _PBPOLRAJAHIOPSPARSE_H
#define _PBPOLRAJAHIOPSPARSE_H

#include <opflow.h>
#include "paramsrajahiop.h"

struct _p_FormPBPOLRAJAHIOPSPARSE {};
typedef struct _p_FormPBPOLRAJAHIOPSPARSE *PBPOLRAJAHIOPSPARSE;

extern PetscErrorCode
OPFLOWSetVariableBoundsArray_PBPOLRAJAHIOPSPARSE(OPFLOW, double *, double *);
extern PetscErrorCode
OPFLOWSetConstraintBoundsArray_PBPOLRAJAHIOPSPARSE(OPFLOW, double *, double *);
extern PetscErrorCode OPFLOWSetInitialGuessArray_PBPOLRAJAHIOPSPARSE(OPFLOW,
                                                                     double *);
extern PetscErrorCode OPFLOWComputeEqualityConstraintsArray_PBPOLRAJAHIOPSPARSE(
    OPFLOW, const double *, double *);
extern PetscErrorCode
OPFLOWComputeInequalityConstraintsArray_PBPOLRAJAHIOPSPARSE(OPFLOW,
                                                            const double *,
                                                            double *);
extern PetscErrorCode
OPFLOWComputeObjectiveArray_PBPOLRAJAHIOPSPARSE(OPFLOW, const double *,
                                                double *);
extern PetscErrorCode
OPFLOWComputeGradientArray_PBPOLRAJAHIOPSPARSE(OPFLOW, const double *,
                                               double *);
extern PetscErrorCode
OPFLOWComputeSparseEqualityConstraintJacobian_PBPOLRAJAHIOPSPARSE(
    OPFLOW, const double *, int *, int *, double *);
extern PetscErrorCode
OPFLOWComputeSparseInequalityConstraintJacobian_PBPOLRAJAHIOPSPARSE(
    OPFLOW, const double *, int *, int *, double *);
extern PetscErrorCode OPFLOWComputeSparseHessian_PBPOLRAJAHIOPSPARSE(
    OPFLOW, const double *, const double *, int *, int *, double *);
extern PetscErrorCode
OPFLOWComputeDenseEqualityConstraintJacobian_PBPOLRAJAHIOPSPARSE(OPFLOW,
                                                                 const double *,
                                                                 double *);
extern PetscErrorCode
OPFLOWComputeDenseInequalityConstraintJacobian_PBPOLRAJAHIOPSPARSE(
    OPFLOW, const double *, double *);
extern PetscErrorCode
OPFLOWComputeDenseHessian_PBPOLRAJAHIOPSPARSE(OPFLOW, const double *,
                                              const double *, double *);

#endif
#endif
#endif
