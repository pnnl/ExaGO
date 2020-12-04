#include <exago_config.h>

/* Embarassingly parallel solver - It does not actually solve the SOPFLOW, but rather
   runs each OPFLOW independently
*/
#ifndef SOPFLOWEMPAR_H
#define SOPFLOWEMPAR_H

typedef struct _p_SOPFLOWSolver_EMPAR *SOPFLOWSolver_EMPAR;

struct _p_SOPFLOWSolver_EMPAR {
  int nxi;
};

#endif
