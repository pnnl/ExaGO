#include <exago_config.h>

#if defined(EXAGO_ENABLE_RAJA)

#ifndef _PBPOLRAJAHIOP_H
#define _PBPOLRAJAHIOP_H

#include <opflow.h>
#include "paramsrajahiop.h"

extern PetscErrorCode
OPFLOWSetVariableBoundsArray_PBPOLRAJAHIOP(OPFLOW, double *, double *);
extern PetscErrorCode
OPFLOWSetConstraintBoundsArray_PBPOLRAJAHIOP(OPFLOW, double *, double *);
extern PetscErrorCode OPFLOWSetInitialGuessArray_PBPOLRAJAHIOP(OPFLOW,
                                                               double *);
extern PetscErrorCode
OPFLOWComputeEqualityConstraintsArray_PBPOLRAJAHIOP(OPFLOW, const double *,
                                                    double *);
extern PetscErrorCode
OPFLOWComputeInequalityConstraintsArray_PBPOLRAJAHIOP(OPFLOW, const double *,
                                                      double *);
extern PetscErrorCode
OPFLOWComputeObjectiveArray_PBPOLRAJAHIOP(OPFLOW, const double *, double *);
extern PetscErrorCode
OPFLOWComputeGradientArray_PBPOLRAJAHIOP(OPFLOW, const double *, double *);
extern PetscErrorCode
OPFLOWComputeSparseEqualityConstraintJacobian_PBPOLRAJAHIOP(OPFLOW,
                                                            const double *,
                                                            int *, int *,
                                                            double *);
extern PetscErrorCode
OPFLOWComputeSparseInequalityConstraintJacobian_PBPOLRAJAHIOP(OPFLOW,
                                                              const double *,
                                                              int *, int *,
                                                              double *);
extern PetscErrorCode
OPFLOWComputeSparseHessian_PBPOLRAJAHIOP(OPFLOW, const double *, const double *,
                                         int *, int *, double *);
extern PetscErrorCode
OPFLOWComputeDenseEqualityConstraintJacobian_PBPOLRAJAHIOP(OPFLOW,
                                                           const double *,
                                                           double *);
extern PetscErrorCode
OPFLOWComputeDenseInequalityConstraintJacobian_PBPOLRAJAHIOP(OPFLOW,
                                                             const double *,
                                                             double *);
extern PetscErrorCode OPFLOWComputeDenseHessian_PBPOLRAJAHIOP(OPFLOW,
                                                              const double *,
                                                              const double *,
                                                              double *);

extern PetscErrorCode
OPFLOWSolutionCallback_PBPOLRAJAHIOP(OPFLOW, const double *, const double *,
                                     const double *, const double *,
                                     const double *, double);
#endif
#endif
