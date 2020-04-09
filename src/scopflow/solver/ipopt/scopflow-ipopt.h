#include <scopflow_config.h>
#if defined(SCOPFLOW_HAVE_IPOPT)

#ifndef SCOPFLOWIPOPT_H
#define SCOPFLOWIPOPT_H

#include <scopflow.h>
#include <IpStdCInterface.h>
#include "../../../opflow/solver/ipopt/aij.h"
#include "../../../opflow/solver/ipopt/sbaij.h"
#include "../../../opflow/solver/ipopt/opflow-ipopt.h"

typedef struct _p_SCOPFLOWSolver_IPOPT *SCOPFLOWSolver_IPOPT;

struct _p_SCOPFLOWSolver_IPOPT {
  
  PetscInt nnz_jac_ge,nnz_jac_gi,nnz_jac_g; /**< Number of nonzeros in the jacobian of the constraints */
  PetscInt nnz_hes; /**< Number of nonzeros in the Lagrangian Hessian */

  Mat      Jac_GeT; /* Transpose of Equality constrained Jacobian matrix */
  Mat      Jac_GiT; /* Transpose of Inequality constrained Jacobian matrix */

  Mat      Hes_sbaij; /* Hessian in symmetric baij format which is needed for the IPOPT solver */
  
  CCMatrix jac_ge;
  CCMatrix jac_gi;
  CCMatrix hes;

  PetscInt *nxi; /* Number of variables for each scenario */
  PetscInt *ngi; /* Number of constraints for each scenario (includes coupling constraints) */
  PetscInt *xstarti; /* Starting location for the variables for scenario i in the big X vector */
  PetscInt *gstarti; /* Starting location for the constraints for scenario i in the big G vector */

  IpoptProblem nlp; /**< Ipopt solver */
  enum ApplicationReturnStatus solve_status;

  PetscScalar obj_factor; /* objective scaling factor set by IPOPT for hessian */
};

#endif
#endif
