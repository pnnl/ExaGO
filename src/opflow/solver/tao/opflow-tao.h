#include <opflow.h>

typedef struct _p_OPFLOWSolver_TAO *OPFLOWSolver_TAO;

struct _p_OPFLOWSolver_TAO {
  
  /* Lagrange multipliers */
  Vec lambda_g;
  Vec lambda_xl;
  Vec lambda_xu;

  Tao nlp; /**< Tao solver */
  PetscBool converged;

  PetscScalar obj_factor; /* objective scaling factor set by IPOPT for hessian */
};
