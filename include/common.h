/*
 * Public header file containing the common objects, API used by different
 * applications
 */

#ifndef COMMON_H
#define COMMON_H

#include <petsc.h>
#include <exago_config.h>

typedef enum { MATPOWER, CSV, JSON, MINIMAL } OutputFormat;

typedef enum { NATIVE = 0, PSSE = 1 } ContingencyFileInputFormat;

typedef enum {
  SOPFLOW_NATIVE_SINGLEPERIOD,
  SOPFLOW_NATIVE_MULTIPERIOD
} ScenarioFileInputFormat;

typedef enum {
  FORECAST_WIND = 1,
  FORECAST_LOAD_P = 2,
  FORECAST_LOAD_Q = 3
} ForecastType;

/**
 * The communicator context
 */
struct _p_COMM {
  MPI_Comm type;    /**< MPI communicator SELF or WORLD */
  PetscMPIInt rank; /**< Process rank */
  PetscMPIInt size; /**< Communicator size */
  PetscInt refct; /**< Reference count to know how many objects are sharing the
                     communicator */
};

typedef struct _p_COMM *COMM;

PETSC_EXTERN PetscErrorCode COMMCreate(MPI_Comm, COMM *);
PETSC_EXTERN PetscErrorCode COMMDestroy(COMM *);

#endif
