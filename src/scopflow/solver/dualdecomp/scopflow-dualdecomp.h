#include <exago_config.h>

/* Embarassingly parallel solver - It does not actually solve the SCOPFLOW, but rather
   runs each OPFLOW independently
*/
#ifndef SCOPFLOWDUALDECOMP_H
#define SCOPFLOWDUALDECOMP_H

typedef struct _p_SCOPFLOWSolver_DUALDECOMP *SCOPFLOWSolver_DUALDECOMP;

struct _p_SCOPFLOWSolver_DUALDECOMP {
  int nxi;
};

#endif
