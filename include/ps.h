/*
 * Public header file containing PS data structs and API
 */

#ifndef __PS_H
#define __PS_H

#include <common.h>
#include <constants.h>
#include <petscdmnetwork.h>

typedef struct _p_PSBUS *PSBUS;
typedef struct _p_PSLOAD *PSLOAD;
typedef struct _p_PSGEN *PSGEN;
typedef struct _p_PSLINE *PSLINE;

/*
 * Type of application that the PS is used for
 */
typedef enum {
  APP_NONE, /**< No application (default) */
  APP_ACPF, /**< AC Power flow */
  APP_ACOPF /**< AC Optimal Power Flow */
} PSApp;

/*
 * Power system public struct
 */
typedef struct _p_PS *PS;

PETSC_EXTERN PetscErrorCode PSReadMatPowerData(PS, const char[]);
PETSC_EXTERN PetscErrorCode PSReadPSSERawData(PS, const char[]);
PETSC_EXTERN PetscErrorCode PSCreate(MPI_Comm, PS *);
PETSC_EXTERN PetscErrorCode PSDestroy(PS *);
PETSC_EXTERN PetscErrorCode PSGENDestroy(PS);
PETSC_EXTERN PetscErrorCode PSSetApplication(PS, void *, PSApp);
PETSC_EXTERN PetscErrorCode PSSetUp(PS);
PETSC_EXTERN PetscErrorCode PSCreateGlobalVector(PS, Vec *);
PETSC_EXTERN PetscErrorCode PSCreateMatrix(PS, Mat *);
PETSC_EXTERN PetscErrorCode PSBUSIsGhosted(PSBUS, PetscBool *);
PETSC_EXTERN PetscErrorCode PSBUSGetNGen(PSBUS, PetscInt *);
PETSC_EXTERN PetscErrorCode PSBUSGetNLoad(PSBUS, PetscInt *);
PETSC_EXTERN PetscErrorCode PSBUSGetGen(PSBUS, PetscInt, PSGEN *);
PETSC_EXTERN PetscErrorCode PSBUSSetGenStatus(PSBUS, char[], PetscInt);
PETSC_EXTERN PetscErrorCode PSBUSGetLoad(PSBUS, PetscInt, PSLOAD *);
PETSC_EXTERN PetscErrorCode PSBUSGetVariableLocation(PSBUS, PetscInt *);
PETSC_EXTERN PetscErrorCode PSBUSGetVariableGlobalLocation(PSBUS, PetscInt *);
PETSC_EXTERN PetscErrorCode PSBUSGetSupportingLines(PSBUS, PetscInt *,
                                                    const PSLINE **);
PETSC_EXTERN PetscErrorCode PSBUSAddShunt(PSBUS, PetscScalar, PetscScalar);
PETSC_EXTERN PetscErrorCode PSLINEGetConnectedBuses(PSLINE, const PSBUS **);
PETSC_EXTERN PetscErrorCode PSLINESetStatus(PSLINE, PetscInt);
PETSC_EXTERN PetscErrorCode PSLINEGetVariableLocation(PSLINE, PetscInt *);
PETSC_EXTERN PetscErrorCode PSLINEGetVariableGlobalLocation(PSLINE, PetscInt *);
PETSC_EXTERN PetscErrorCode PSGetNumActiveLines(PS, PetscInt *, PetscInt *);
PETSC_EXTERN PetscErrorCode PSGetNumGenerators(PS, PetscInt *, PetscInt *);
PETSC_EXTERN PetscErrorCode PSGetNumLoads(PS, PetscInt *, PetscInt *);
PETSC_EXTERN PetscErrorCode PSGetNumActiveGenerators(PS, PetscInt *,
                                                     PetscInt *);
PETSC_EXTERN PetscErrorCode PSGetNumGlobalBuses(PS, PetscInt *);
PETSC_EXTERN PetscErrorCode PSGENSetStatus(PSGEN, PetscInt);
PETSC_EXTERN PetscErrorCode PSLOADSetStatus(PSLOAD, PetscInt);
PETSC_EXTERN PetscErrorCode PSCheckandSetRefBus(PS);
PETSC_EXTERN PetscErrorCode PSGetTotalGeneration(PS, PetscScalar *,
                                                 PetscScalar *);
PETSC_EXTERN PetscErrorCode PSComputeParticipationFactors(PS);
PETSC_EXTERN PetscErrorCode PSSetGenStatus(PS, PetscInt, const char *,
                                           PetscInt);
PETSC_EXTERN PetscErrorCode PSSetLineStatus(PS, PetscInt, PetscInt,
                                            const char *, PetscInt);
PETSC_EXTERN PetscErrorCode PSCopy(PS, PS);

/*
  PSSetGenPowerLimits - Sets generator real and reactive power limits given a
bus number and generator id

  Input Parameters
+ ps   - the PS object
. gbus - bus number for the generator
. id   - generator id
. pt   - Max. generator real power capacity (MW)
. pb   - Min. generator real power capacity (MW)
. qt   - Max. generator reactive power capacity (MVAr)
- qb   - Min. generator reactive power capacity (MVAr)

  Notes: Use EXAGO_IGNORE to not set the value
*/
PETSC_EXTERN PetscErrorCode PSSetGenPowerLimits(PS, PetscInt, const char *,
                                                PetscScalar, PetscScalar,
                                                PetscScalar, PetscScalar);

/*
  PSGetGenPowerLimits - Returns generator real and reactive power limits given a
bus number and generator id

  Input Parameters
+ ps   - the PS object
. gbus - bus number for the generator
- id   - generator id

  Output Paramters
+ pt   - Max. generator real power capacity
. pb   - Min. generator real power capacity
. qt   - Max. generator reactive power capacity
- qb   - Min. generator reactive power capacity

  Notes: Use NULL to ignore retrieving value
*/
PETSC_EXTERN PetscErrorCode PSGetGenPowerLimits(PS, PetscInt, const char *,
                                                PetscScalar *, PetscScalar *,
                                                PetscScalar *, PetscScalar *);

/*
  PSSetGenStatus - Sets generator status given the bus name and id
  Input Parameters
+ ps     - the PS object
. gbus   - generator bus
. gid    - generator id
- status - Generator status (1 = ON, 0 = OFF)

  Notes: PSSetUp() must be called before calling PSSetGenStatus()
*/
PETSC_EXTERN PetscErrorCode PSSetGenStatus(PS, PetscInt, const char *,
                                           PetscInt);

/*
  PSSetGenStatus - Gets generator status given the bus name and id

  Input Parameters
+ ps     - the PS object
. gbus   - generator bus
- gid    - generator id

  Output Parameters
. status - Generator status

  Notes: PSSetUp() must be called before calling PSSetGenStatus()
*/
PETSC_EXTERN PetscErrorCode PSGetGenStatus(PS, PetscInt, const char *,
                                           PetscInt *);

PETSC_EXTERN PetscErrorCode PSSetLineStatus(PS, PetscInt, PetscInt,
                                            const char *, PetscInt);

/*
  PSGetGenDispatch - Returns generator real (MW) and reactive power (MVAR)
dispatch given a bus number and generator id

  Input Parameters
+ ps   - the PS object
. gbus - bus number for the generator
- id   - generator id

  Output Paramters
+ pg   - Real power dispatch (MW)
- qg   - Reactive power dispatch (MVAr)

  Notes: This function should be called after OPFLOWSolutionToPS()
*/
PETSC_EXTERN PetscErrorCode PSGetGenDispatch(PS, PetscInt, const char *,
                                             PetscScalar *, PetscScalar *);

#endif
