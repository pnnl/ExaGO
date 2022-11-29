/*
 * Public header file for optimal power flow application.
 */

#ifndef OPFLOW_H
#define OPFLOW_H

#include <ps.h>
#include <utils.h>

/* Models */
#define OPFLOWMODEL_PBPOL "POWER_BALANCE_POLAR"
#define OPFLOWMODEL_PBCAR "POWER_BALANCE_CARTESIAN"
#define OPFLOWMODEL_IBCAR "CURRENT_BALANCE_CARTESIAN"
#define OPFLOWMODEL_IBCAR2 "CURRENT_BALANCE_CARTESIAN2"
#define OPFLOWMODEL_PBPOLHIOP "POWER_BALANCE_HIOP"
#define OPFLOWMODEL_PBPOLRAJAHIOP "PBPOLRAJAHIOP"
#define OPFLOWMODEL_DCOPF "DCOPF"

/* Solvers */
#define OPFLOWSOLVER_IPOPT "IPOPT"
#define OPFLOWSOLVER_HIOP "HIOP"
#define OPFLOWSOLVER_HIOPSPARSE "HIOPSPARSE"

typedef struct _p_OPFLOW *OPFLOW;

typedef enum {
#define OPFLOW_OBJ_DEF(x) x,
#include "private/opflow_enum.def"
#undef OPFLOW_OBJ_DEF
} OPFLOWObjectiveType;

/* Generator bus voltage type */
typedef enum {
#define OPFLOW_GBV_DEF(x) x,
#include "private/opflow_enum.def"
#undef OPFLOW_GBV_DEF
} OPFLOWGenBusVoltageType;

typedef enum {
#define OPFLOW_INIT_DEF(x) x,
#include "private/opflow_enum.def"
#undef OPFLOW_INIT_DEF
} OPFLOWInitializationType;

/**
 * \brief Option declarations and default values.
 *
 * The values are printed in the help messages here:
 * \see src/utils/utils.cpp
 *
 * The default values are used and the options are declared via PETSc here:
 * \see src/opflow/interface/opflow.cpp
 */
namespace OPFLOWOptions {

/*
 * Default model/solver combination is a valid Ipopt configuration, and falls
 * back to HiOp. If niether is enabled, error at build-time (although the CMake
 * should error at configure-time anyways).
 */
#if defined(EXAGO_ENABLE_IPOPT)
const auto model = ExaGOStringOption("-opflow_model", "OPFLOW model name",
                                     "POWER_BALANCE_POLAR",
                                     {
#ifdef EXAGO_ENABLE_HIOP
                                         "POWER_BALANCE_HIOP",
                                         "PBPOLRAJAHIOP",
#endif
                                     });
const auto solver =
    ExaGOStringOption("-opflow_solver", "OPFLOW solver type", "IPOPT",
                      {
#ifdef EXAGO_ENABLE_HIOP
                          "HIOP",
                          "HIOPSPARSE",
#endif
                      });
#elif defined(EXAGO_ENABLE_HIOP)
const auto model = ExaGOStringOption("-opflow_model", "OPFLOW model name",
                                     "POWER_BALANCE_HIOP",
                                     {
                                         "POWER_BALANCE_POLAR",
                                         "PBPOLRAJAHIOP",
                                     });
const auto solver = ExaGOStringOption("-opflow_solver", "OPFLOW solver type",
                                      "OPFLOWSOLVER_HIOP",
                                      {
#ifdef EXAGO_ENABLE_IPOPT
                                          "IPOPT",
#endif
                                          "HIOPSPARSE",
                                      });
#else
#error "ExaGO must be built with at least one solver enabled"
#endif

const auto initialization =
    ExaGOStringOption("-opflow_initialization", "Type of OPFLOW initialization",
                      "OPFLOWINIT_MIDPOINT",
                      {
                          "OPFLOWINIT_FROMFILE",
                          "OPFLOWINIT_ACPF",
                          "OPFLOWINIT_FLATSTART",
                          "OPFLOWINIT_DCOPF",
                      });
const auto objective = ExaGOStringOption(
    "-opflow_objective", "Type of OPFLOW objective", "MIN_GEN_COST",
    {
        "MIN_GENSETPOINT_DEVIATION",
        "NO_OBJ",
    });
const auto genbusvoltage = ExaGOStringOption(
    "-opflow_genbusvoltage", "Type of OPFLOW gen bus voltage control",
    "VARIABLE_WITHIN_BOUNDS",
    {
        "FIXED_WITHIN_QBOUNDS",
        "FIXED_AT_SETPOINT",
    });
const auto has_gensetpoint =
    ExaGOBoolOption("-opflow_has_gensetpoint",
                    "Use set-points for generator real power", PETSC_FALSE);
const auto use_agc = ExaGOBoolOption(
    "-opflow_use_agc", "Use automatic generation control (AGC)", PETSC_FALSE);
const auto tolerance =
    ExaGORealOption("-opflow_tolerance", "Optimization tolerance", 1e-6);

const auto ignore_lineflow_constraints =
    ExaGOBoolOption("-opflow_ignore_lineflow_constraints",
                    "Ignore line flow constraints?", PETSC_FALSE);
const auto include_loadloss_variables =
    ExaGOBoolOption("-opflow_include_loadloss_variables",
                    "Ignore line flow constraints?", PETSC_FALSE);
const auto loadloss_penalty =
    ExaGORealOption("-opflow_loadloss_penalty", "Penalty for load loss", 1e3);
const auto include_powerimbalance_variables =
    ExaGOBoolOption("-opflow_include_powerimbalance_variables",
                    "Allow power imbalance?", PETSC_FALSE);
const auto powerimbalance_penalty = ExaGORealOption(
    "-opflow_powerimbalance_penalty", "Power imbalance penalty", 1e4);

#ifdef EXAGO_ENABLE_HIOP
const auto hiop_compute_mode =
    ExaGOStringOption("-hiop_compute_mode", "Set compute mode for HiOp solver",
                      "auto", {"cpu", "hybrid", "gpu"});
const auto hiop_verbosity_level =
    ExaGOIntOption("-hiop_verbosity_level",
                   "Set verbosity level for HiOp solver, between 0 and 12", 0);
const auto hiop_mem_space =
    ExaGOStringOption("-hiop_mem_space", "Set memory space for HiOp solver",
                      "host", {"default", "host", "um", "device"});

#ifdef EXAGO_ENABLE_IPOPT
const auto hiop_ipopt_debug = ExaGOBoolOption(
    "-hiop_ipopt_debug", "Flag enabling debugging HIOP code with IPOPT",
    PETSC_FALSE);
#endif

#endif

} // namespace OPFLOWOptions

PETSC_EXTERN PetscErrorCode OPFLOWSetModel(OPFLOW, const char[]);
PETSC_EXTERN PetscErrorCode OPFLOWGetModel(OPFLOW, char *);
PETSC_EXTERN PetscErrorCode OPFLOWSetSolver(OPFLOW, const char[]);
PETSC_EXTERN PetscErrorCode OPFLOWGetSolver(OPFLOW, char *);
PETSC_EXTERN PetscErrorCode
OPFLOWModelRegister(OPFLOW, const char[], PetscErrorCode (*create)(OPFLOW));
PETSC_EXTERN PetscErrorCode
OPFLOWSolverRegister(OPFLOW, const char[], PetscErrorCode (*create)(OPFLOW));
PETSC_EXTERN PetscErrorCode OPFLOWCreate(MPI_Comm, OPFLOW *);
PETSC_EXTERN PetscErrorCode OPFLOWDestroy(OPFLOW *);
PETSC_EXTERN PetscErrorCode OPFLOWReadMatPowerData(OPFLOW, const char[]);
PETSC_EXTERN PetscErrorCode OPFLOWSetUp(OPFLOW);
PETSC_EXTERN PetscErrorCode OPFLOWCreateGlobalVector(OPFLOW, Vec *);
PETSC_EXTERN PetscErrorCode OPFLOWCreateMatrix(OPFLOW, Mat *);
PETSC_EXTERN PetscErrorCode OPFLOWSolve(OPFLOW);
PETSC_EXTERN PetscErrorCode OPFLOWSetInitialGuess(OPFLOW, Vec, Vec);
PETSC_EXTERN PetscErrorCode OPFLOWSetTolerance(OPFLOW, PetscReal);
PETSC_EXTERN PetscErrorCode OPFLOWSetHIOPComputeMode(OPFLOW, const char *);
PETSC_EXTERN PetscErrorCode OPFLOWGetHIOPComputeMode(OPFLOW, char *);
PETSC_EXTERN PetscErrorCode OPFLOWSetHIOPMemSpace(OPFLOW, const char *);
PETSC_EXTERN PetscErrorCode OPFLOWGetHIOPMemSpace(OPFLOW, char *);
PETSC_EXTERN PetscErrorCode OPFLOWSetHIOPVerbosityLevel(OPFLOW, int);
PETSC_EXTERN PetscErrorCode OPFLOWGetHIOPVerbosityLevel(OPFLOW, int *);
PETSC_EXTERN PetscErrorCode OPFLOWGetTolerance(OPFLOW, PetscReal *);
PETSC_EXTERN PetscErrorCode OPFLOWGetNumIterations(OPFLOW, PetscInt *);
PETSC_EXTERN PetscErrorCode OPFLOWGetObjectiveType(OPFLOW,
                                                   OPFLOWObjectiveType *);
PETSC_EXTERN PetscErrorCode OPFLOWSetObjectiveType(OPFLOW, OPFLOWObjectiveType);
PETSC_EXTERN PetscErrorCode OPFLOWGetObjective(OPFLOW, PetscReal *);
PETSC_EXTERN PetscErrorCode OPFLOWGetVariableBounds(OPFLOW, Vec *, Vec *);
PETSC_EXTERN PetscErrorCode OPFLOWGetConstraintJacobian(OPFLOW, Mat *, Mat *);
PETSC_EXTERN PetscErrorCode OPFLOWGetHessian(OPFLOW, Mat *, PetscScalar *);
PETSC_EXTERN PetscErrorCode OPFLOWGetSolution(OPFLOW, Vec *);
PETSC_EXTERN PetscErrorCode OPFLOWGetConvergenceStatus(OPFLOW, PetscBool *);
PETSC_EXTERN PetscErrorCode OPFLOWGetConstraints(OPFLOW, Vec *);
PETSC_EXTERN PetscErrorCode OPFLOWGetConstraintMultipliers(OPFLOW, Vec *);

PETSC_EXTERN PetscErrorCode OPFLOWComputeVariableBounds(OPFLOW, Vec, Vec);
PETSC_EXTERN PetscErrorCode OPFLOWComputeGradient(OPFLOW, Vec, Vec);
PETSC_EXTERN PetscErrorCode OPFLOWComputeGradientArray(OPFLOW, const double *,
                                                       double *);
PETSC_EXTERN PetscErrorCode OPFLOWComputeObjective(OPFLOW, Vec, PetscReal *);
PETSC_EXTERN PetscErrorCode OPFLOWComputeObjectiveArray(OPFLOW, const double *,
                                                        double *);
PETSC_EXTERN PetscErrorCode OPFLOWComputeConstraints(OPFLOW, Vec, Vec);
PETSC_EXTERN PetscErrorCode OPFLOWComputeEqualityConstraints(OPFLOW, Vec, Vec);
PETSC_EXTERN PetscErrorCode OPFLOWComputeInequalityConstraints(OPFLOW, Vec,
                                                               Vec);
PETSC_EXTERN PetscErrorCode
OPFLOWComputeEqualityConstraintsArray(OPFLOW, const double *, double *);
PETSC_EXTERN PetscErrorCode
OPFLOWComputeInequalityConstraintsArray(OPFLOW, const double *, double *);

PETSC_EXTERN PetscErrorCode OPFLOWGetConstraintBounds(OPFLOW, Vec *, Vec *);
PETSC_EXTERN PetscErrorCode OPFLOWComputeConstraintBounds(OPFLOW, Vec, Vec);
PETSC_EXTERN PetscErrorCode OPFLOWComputeConstraintJacobian(OPFLOW, Vec, Mat,
                                                            Mat);
PETSC_EXTERN PetscErrorCode OPFLOWComputeHessian(OPFLOW, Vec, Vec, PetscScalar,
                                                 Mat);

PETSC_EXTERN PetscErrorCode OPFLOWPrintSolution(OPFLOW);
PETSC_EXTERN PetscErrorCode OPFLOWSaveSolution(OPFLOW, OutputFormat,
                                               const char[]);

PETSC_EXTERN PetscErrorCode OPFLOWGetVariableOrdering(OPFLOW, int **);
PETSC_EXTERN PetscErrorCode OPFLOWGetSizes(OPFLOW, int *, int *, int *);

PETSC_EXTERN PetscErrorCode OPFLOWHasGenSetPoint(OPFLOW, PetscBool);
PETSC_EXTERN PetscErrorCode OPFLOWGetHasGenSetPoint(OPFLOW, PetscBool *);
PETSC_EXTERN PetscErrorCode OPFLOWHasLoadLoss(OPFLOW, PetscBool);
PETSC_EXTERN PetscErrorCode OPFLOWGetHasLoadLoss(OPFLOW, PetscBool *);
PETSC_EXTERN PetscErrorCode OPFLOWHasBusPowerImbalance(OPFLOW, PetscBool);
PETSC_EXTERN PetscErrorCode OPFLOWGetHasBusPowerImbalance(OPFLOW, PetscBool *);
PETSC_EXTERN PetscErrorCode OPFLOWUseAGC(OPFLOW, PetscBool);
PETSC_EXTERN PetscErrorCode OPFLOWGetUseAGC(OPFLOW, PetscBool *);

PETSC_EXTERN PetscErrorCode OPFLOWSetGenBusVoltageType(OPFLOW,
                                                       OPFLOWGenBusVoltageType);
PETSC_EXTERN PetscErrorCode
OPFLOWGetGenBusVoltageType(OPFLOW, OPFLOWGenBusVoltageType *);
PETSC_EXTERN PetscErrorCode
    OPFLOWSetInitializationType(OPFLOW, OPFLOWInitializationType);
PETSC_EXTERN PetscErrorCode
OPFLOWGetInitializationType(OPFLOW, OPFLOWInitializationType *);
PETSC_EXTERN PetscErrorCode OPFLOWIgnoreLineflowConstraints(OPFLOW, PetscBool);
PETSC_EXTERN PetscErrorCode OPFLOWGetIgnoreLineflowConstraints(OPFLOW,
                                                               PetscBool *);
PETSC_EXTERN PetscErrorCode OPFLOWSetLoadLossPenalty(OPFLOW, PetscReal);
PETSC_EXTERN PetscErrorCode OPFLOWGetLoadLossPenalty(OPFLOW, PetscReal *);
PETSC_EXTERN PetscErrorCode OPFLOWSetBusPowerImbalancePenalty(OPFLOW,
                                                              PetscReal);
PETSC_EXTERN PetscErrorCode OPFLOWGetBusPowerImbalancePenalty(OPFLOW,
                                                              PetscReal *);
PETSC_EXTERN PetscErrorCode OPFLOWSetWeight(OPFLOW, PetscScalar);

PETSC_EXTERN PetscErrorCode OPFLOWSkipOptions(OPFLOW, PetscBool);

/* OPFLOWGetPS - Gets the underlying PS object

  Input Parameters
. opflow - the OPFLOW object

  Output Parameters
. ps - the ps object

  Notes: This function returns the PS object that holds the network data. Using
the PS object one can make changes to the network parameters. A typical case is
         changing some network parameters before solving opflow.
         OPFLOWSetUpPS() must be called before OPFLOWGetPS()
*/
PETSC_EXTERN PetscErrorCode OPFLOWGetPS(OPFLOW, PS *);

/* OPFLOWSetUpPS - Sets the underlying PS network object to be used by OPFLOW

  Input Parameters
. opflow - the OPFLOW object

  Notes: This function is an intermediate function that can be called for
setting up the PS network object prior to solving OPFLOW. A typical use-case is
some network parameter needs changing before solving opflow. In such case, the
work flow would be
  1. OPFLOWCreate();
  2. OPFLOWReadMatPowerData();
  3. OPFLOWSetUpPS();
  4. OPFLOWGetPS();
  ... change the network data by the PS object retrieved via OPFLOWGetPS().
  5. OPFLOWSolve();
 Skip steps 3 and 4 if no network changes are needed to be done.
*/
PETSC_EXTERN PetscErrorCode OPFLOWSetUpPS(OPFLOW);

/*
  OPFLOWSolutionToPS - Updates the PS struct from OPFLOW solution

  Input Parameters:
. opflow - the OPFLOW object

  Notes: Updates the different fields in the PS struct from the OPFLOW solution.
         This function should be called after OPFLOWSolve() before retrieving
         values from the PS struct.
*/
PETSC_EXTERN PetscErrorCode OPFLOWSolutionToPS(OPFLOW);

/*
  OPFLOWSetLinesMonitored - List of lines to monitor. The flows for these lines
                            are included as inequality constraints in OPFLOW

 Input Parameter:
+ opflow      - OPFLOW object
. nkvlevels   - Number of kvlevels to monitor (Use -1 to monitor all kvlevels)
. kvlevels    - line kvlevels to monitor
- monitorfile - File with list of lines to monitor.

  Notes:
    The lines to monitor are either specified through a file OR by
    kvlevels, but not both. Use NULL for monitorfile if file is not set.
    If monitorfile is given then the kvlevels are ignored.

    Lines are specified in the file in the format frombus,tobus where
    frombus and tobus are the from and to bus numbers for the line.

    This function should be called after OPFLOWSetupPS() is called
*/
PETSC_EXTERN PetscErrorCode OPFLOWSetLinesMonitored(OPFLOW, PetscInt,
                                                    const PetscScalar *,
                                                    const char *);

typedef PetscErrorCode (*OPFLOWAuxObjectiveFunction)(OPFLOW, const double *,
                                                     double *, void *);
typedef PetscErrorCode (*OPFLOWAuxGradientFunction)(OPFLOW, const double *,
                                                    double *, void *);
typedef PetscErrorCode (*OPFLOWAuxHessianFunction)(OPFLOW, const double *, Mat,
                                                   void *);

PETSC_EXTERN PetscErrorCode OPFLOWSetAuxillaryObjective(
    OPFLOW, OPFLOWAuxObjectiveFunction, OPFLOWAuxGradientFunction,
    OPFLOWAuxHessianFunction, void *);

PETSC_EXTERN PetscErrorCode OPFLOWSetUpdateVariableBoundsFunction(
    OPFLOW, PetscErrorCode (*)(OPFLOW, Vec, Vec, void *), void *);

#endif
