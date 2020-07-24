#if defined(EXAGO_HAVE_PIPS)

#ifndef SCOPFLOWPIPS_H
#define SCOPFLOWPIPS_H

#include <scopflow.h>
#include <Drivers/parallelPipsNlp_C_Callback.h>
#include <../src/mat/impls/aij/seq/aij.h>
#include <../src/mat/impls/sbaij/seq/sbaij.h>
#include "../../../opflow/solver/ipopt/opflow-ipopt.h"

typedef struct _p_SCOPFLOWSolver_PIPS *SCOPFLOWSolver_PIPS;

struct _p_SCOPFLOWSolver_PIPS {
  
  PetscInt nnz_jac_ge,nnz_jac_gi,nnz_jac_g; /**< Number of nonzeros in the jacobian of the constraints */
  PetscInt nnz_hes; /**< Number of nonzeros in the Lagrangian Hessian */

  PetscInt *nxi; /* Number of variables for each scenario */
  PetscInt *ngi; /* Number of constraints for each scenario (includes coupling constraints) */
  PetscInt *xstarti; /* Starting location for the variables for scenario i in the big X vector */
  PetscInt *gstarti; /* Starting location for the constraints for scenario i in the big G vector */

  Mat      *Jac_GicoupT; /* Transpose of coupling jacobian */
  PipsNlpProblemStructPtr nlp;

};

#endif
#endif
