#include <exago_config.h>
#if defined(EXAGO_HAVE_IPOPT)

#ifndef SCOPFLOWIPOPT_H
#define SCOPFLOWIPOPT_H

#include <scopflow.h>
#include <IpStdCInterface.h>
#include "../../../opflow/solver/ipopt/opflow-ipopt.h"

typedef struct _p_SCOPFLOWSolver_IPOPT *SCOPFLOWSolver_IPOPT;

struct _p_SCOPFLOWSolver_IPOPT {
  
  PetscInt nnz_jac_ge,nnz_jac_gi,nnz_jac_g; /**< Number of nonzeros in the jacobian of the constraints */
  PetscInt nnz_hes; /**< Number of nonzeros in the Lagrangian Hessian */

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
