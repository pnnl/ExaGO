#include <exago_config.h>
#if defined(EXAGO_ENABLE_IPOPT)

#ifndef SOPFLOWIPOPT_H
#define SOPFLOWIPOPT_H

#include <sopflow.h>
#include <IpStdCInterface.h>
#include "../../../opflow/solver/ipopt/opflow-ipopt.h"

typedef struct _p_SOPFLOWSolver_IPOPT *SOPFLOWSolver_IPOPT;

struct _p_SOPFLOWSolver_IPOPT {
  
  PetscInt nnz_jac_g; /**< Number of nonzeros in the jacobian of the constraints */
  PetscInt nnz_hes; /**< Number of nonzeros in the Lagrangian Hessian */

  IpoptProblem nlp; /**< Ipopt solver */
  enum ApplicationReturnStatus solve_status;

  PetscScalar obj_factor; /* objective scaling factor set by IPOPT for hessian */
};

#endif
#endif
