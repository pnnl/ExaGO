#if defined(SCOPFLOW_HAVE_IPOPT)

#include <opflow.h>
#include <IpStdCInterface.h>

typedef struct _p_OPFLOWSolver_IPOPT *OPFLOWSolver_IPOPT;

struct _p_OPFLOWSolver_IPOPT {
  
  PetscInt nnz_jac_g; /**< Number of nonzeros in the jacobian of the constraints */
  PetscInt nnz_hes; /**< Number of nonzeros in the Lagrangian Hessian */

  /* Lagrange multipliers */
  Vec lambda_g;
  Vec lambda_xl;
  Vec lambda_xu;

  IpoptProblem nlp; /**< Ipopt solver */
  enum ApplicationReturnStatus solve_status;

  PetscScalar obj_factor; /* objective scaling factor set by IPOPT for hessian */
};

#endif
