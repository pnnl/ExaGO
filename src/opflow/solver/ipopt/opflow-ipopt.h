#if defined(SCOPFLOW_HAVE_IPOPT)

#include <opflow.h>
#include <IpStdCInterface.h>

/* Data structure for converting matrix in PETSc format (compressed row aij) to matrix market format (row idx, col idx, value)
*/
struct _p_CCMatrix {
  PetscInt    *rowidx;
  PetscInt    *colptr;
  PetscScalar *values;
};

typedef struct _p_CCMatrix *CCMatrix;

typedef struct _p_OPFLOWSolver_IPOPT *OPFLOWSolver_IPOPT;

struct _p_OPFLOWSolver_IPOPT {
  
  PetscInt nnz_jac_g; /**< Number of nonzeros in the jacobian of the constraints */
  PetscInt nnz_hes; /**< Number of nonzeros in the Lagrangian Hessian */

  Mat      Jac_GeT; /* Transpose of Equality constrained Jacobian matrix */
  Mat      Jac_GiT; /* Transpose of Inequality constrained Jacobian matrix */
  
  CCMatrix jac_ge;
  CCMatrix jac_gi;

  /* Lagrange multipliers */
  Vec lambda_g;
  Vec lambda_xl;
  Vec lambda_xu;

  IpoptProblem nlp; /**< Ipopt solver */
  enum ApplicationReturnStatus solve_status;

  PetscScalar obj_factor; /* objective scaling factor set by IPOPT for hessian */
};

#endif
