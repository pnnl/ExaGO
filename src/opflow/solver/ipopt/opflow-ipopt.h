#include <exago_config.h>
#if defined(EXAGO_ENABLE_IPOPT)

#ifndef OPFLOWIPOPT_H
#define OPFLOWIPOPT_H

#include <opflow.h>
#include <IpStdCInterface.h>

typedef struct _p_OPFLOWSolver_IPOPT *OPFLOWSolver_IPOPT;

struct _p_OPFLOWSolver_IPOPT {
  
  PetscInt nnz_jac_ge,nnz_jac_gi,nnz_jac_g; /**< Number of nonzeros in the jacobian of the constraints */
  PetscInt nnz_hes; /**< Number of nonzeros in the Lagrangian Hessian */

  IpoptProblem nlp; /**< Ipopt solver */
  enum ApplicationReturnStatus solve_status;

  PetscScalar obj_factor; /* objective scaling factor set by IPOPT for hessian */
};
#endif
#endif
