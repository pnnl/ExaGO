#include <scopflow_config.h>

/* Embarassingly parallel solver - It does not actually solve the SCOPFLOW, but rather
   runs each OPFLOW independently
*/
#ifndef SCOPFLOWEMPAR_H
#define SCOPFLOWEMPAR_H

typedef struct _p_SCOPFLOWSolver_EMPAR *SCOPFLOWSolver_EMPAR;

struct _p_SCOPFLOWSolver_EMPAR {
  int nxi;
};

#endif
