#include <opflow.h>

typedef struct _p_OPFLOWSolver_TAO *OPFLOWSolver_TAO;

struct _p_OPFLOWSolver_TAO {

  Tao nlp; /**< Tao solver */
  PetscBool converged;

  Vec Glineq, Guineq; /* lower and upper bounds on inequality constraints */
};
