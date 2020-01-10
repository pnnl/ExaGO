#include <opflow.h>

typedef struct _p_OPFLOWSolver_TAO *OPFLOWSolver_TAO;

struct _p_OPFLOWSolver_TAO {
  
  Tao nlp; /**< Tao solver */
  PetscBool converged;

  PetscScalar obj_factor; /* objective scaling factor set by IPOPT for hessian */
};
