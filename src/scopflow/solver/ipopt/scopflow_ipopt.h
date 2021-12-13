#include <exago_config.h>
#if defined(EXAGO_ENABLE_IPOPT)

#ifndef SCOPFLOWIPOPT_H
#define SCOPFLOWIPOPT_H

#include "../../../opflow/solver/ipopt/opflow_ipopt.h"
#include <IpStdCInterface.h>
#include <scopflow.h>

typedef struct _p_SCOPFLOWSolver_IPOPT *SCOPFLOWSolver_IPOPT;

struct _p_SCOPFLOWSolver_IPOPT {

  PetscInt
      nnz_jac_g; /**< Number of nonzeros in the jacobian of the constraints */
  PetscInt nnz_hes; /**< Number of nonzeros in the Lagrangian Hessian */

  IpoptProblem nlp; /**< Ipopt solver */
  enum ApplicationReturnStatus solve_status;

  PetscScalar
      obj_factor; /* objective scaling factor set by IPOPT for hessian */
};

#endif
#endif
