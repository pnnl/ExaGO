/*
 * Public header file for security constrained optimal power flow application.
 *
 */

#ifndef SCOPFLOW_H
#define SCOPFLOW_H

#include <opflow.h>
#include <ps.h>

/* Models */
#define SCOPFLOWMODEL_GENRAMP "GENRAMP"
/* Model for multi-period SCOPFLOW */
#define SCOPFLOWMODEL_GENRAMPT "GENRAMPT"

/* Solvers */
#define SCOPFLOWSOLVER_IPOPT "IPOPT"
/* Embarassingly parallel solver - solves each OPFLOW independently */
#define SCOPFLOWSOLVER_EMPAR "EMPAR"
/* Primal decomposition-based HiOp solver */
#define SCOPFLOWSOLVER_HIOP "HIOP"

/* Initialization and Parameters*/
#define SCOPFLOW_INITIALIZATION "ACPF"
#define SCOPFLOW_GENBUSVOLTAGE "VARIABLE_WITHIN_BOUNDS"

typedef struct _p_SCOPFLOW *SCOPFLOW;

namespace SCOPFLOWOptions {

const auto model = ExaGOStringOption("-scopflow_model", "SCOPFLOW model type",
                                     "GENRAMP", {"GENRAMPT"});

const auto solver =
    ExaGOStringOption("-scopflow_solver", "SCOPFLOW solver type",
#if defined(EXAGO_ENABLE_IPOPT) && defined(EXAGO_ENABLE_HIOP)
                      "IPOPT", {"EMPAR", "HIOP"});
#else
#if defined(EXAGO_ENABLE_IPOPT)
                      "IPOPT", {"EMPAR"});
#elif defined(EXAGO_ENABLE_HIOP)
                      "HIOP", {"EMPAR"});
#else
#error "At least one solver must be enabled!"
#endif
#endif

const auto subproblem_model = ExaGOStringOption(
    "-scopflow_subproblem_model", "SCOPFLOW subproblem model type",
    OPFLOWOptions::model.default_value, OPFLOWOptions::model.possible_values);
const auto subproblem_solver = ExaGOStringOption(
    "-scopflow_subproblem_solver", "SCOPFLOW subproblem solver type",
    OPFLOWOptions::solver.default_value, OPFLOWOptions::solver.possible_values);
#ifdef EXAGO_ENABLE_HIOP
const auto compute_mode =
    ExaGOStringOption("-hiop_compute_mode", "SCOPFLOW subproblem compute mode",
                      OPFLOWOptions::hiop_compute_mode.default_value,
                      OPFLOWOptions::hiop_compute_mode.possible_values);
const auto verbosity_level = ExaGOIntOption(
    "-hiop_verbosity_level", "SCOPFLOW subproblem verbosity level",
    OPFLOWOptions::hiop_verbosity_level.default_value);
#endif
const auto iscoupling = ExaGOBoolOption(
    "-scopflow_iscoupling",
    "Include coupling between first stage and second stage", PETSC_TRUE);

const auto Nc =
    ExaGOIntOption("-scopflow_Nc", "Number of second-stage scenarios", 0);
const auto mode = ExaGOIntOption(
    "-scopflow_mode", "Operation mode: Preventive (0) or Corrective (1)", 1);
const auto enable_multiperiod = ExaGOBoolOption(
    "-scopflow_enable_multiperiod", "Multi-period SCOPFLOW?", PETSC_FALSE);
const auto enable_powerimbalance =
    ExaGOBoolOption("-opflow_include_powerimbalance_variables",
                    "Allow power imbalance", PETSC_FALSE);
const auto ignore_lineflow_constraints = ExaGOBoolOption(
    "-opflow_ignore_lineflow_constraints", "Allow power imbalance",
    OPFLOWOptions::ignore_lineflow_constraints.default_value);

const auto tolerance =
    ExaGORealOption("-scopflow_tolerance", "Optimization tolerance",
                    OPFLOWOptions::tolerance.default_value);

namespace {
static const std::string onlymultiperiod =
    " (only when multiperiod is enabled)";
}
/* Options only relevant when multiperiod is enabled */
const auto windgenprofile = ExaGOStringOption(
    "-scopflow_windgenprofile", "Wind generation file" + onlymultiperiod,
    "/path/to/windgenprofile", {});
const auto ploadprofile = ExaGOStringOption(
    "-scopflow_ploadprofile", "Active power load profile" + onlymultiperiod,
    "/path/to/ploadprofile", {});
const auto qloadprofile = ExaGOStringOption(
    "-scopflow_qloadprofile", "Reactive power load profile" + onlymultiperiod,
    "/path/to/qloadprofile", {});

const auto dT = ExaGORealOption(
    "-scopflow_dT", "Length of time-step (minutes)" + onlymultiperiod, 1.0);
const auto duration = ExaGORealOption(
    "-scopflow_duration", "Time horizon (hours)" + onlymultiperiod, 1.0);

} // namespace SCOPFLOWOptions

PETSC_EXTERN PetscErrorCode SCOPFLOWEnableMultiPeriod(SCOPFLOW, PetscBool);
PETSC_EXTERN PetscErrorCode SCOPFLOWSetModel(SCOPFLOW, const char[]);
PETSC_EXTERN PetscErrorCode SCOPFLOWModelRegister(
    SCOPFLOW, const char[], PetscErrorCode (*create)(SCOPFLOW));
PETSC_EXTERN PetscErrorCode SCOPFLOWSetSolver(SCOPFLOW, const char[]);
PETSC_EXTERN PetscErrorCode SCOPFLOWCreate(MPI_Comm, SCOPFLOW *);
PETSC_EXTERN PetscErrorCode SCOPFLOWDestroy(SCOPFLOW *);
PETSC_EXTERN PetscErrorCode SCOPFLOWSetNetworkData(SCOPFLOW, const char[]);

PETSC_EXTERN PetscErrorCode
SCOPFLOWSetInitilizationType(SCOPFLOW, OPFLOWInitializationType type);
PETSC_EXTERN PetscErrorCode
    SCOPFLOWSetGenBusVoltageType(SCOPFLOW, OPFLOWGenBusVoltageType);

/*
 * The native format for Contingency input file is as follows:
 * Notes: Each field in the contingency file has the following format
 * Num,Type,Bus,Fbus,Tbus,Id,Status,prob
 * Num - Contingency number
 * Type - Type of contingency (Generator, Branch, Transformer, Load)
 * Bus - The equipment bus number
 * Fbus - From bus number (only for branch and transformer contingencies)
 * Tbus - To bus number (only for branch and transformer contingencies)
 * Id   - The equipment ID (2-char string)
 * Status - The status to be set for the outaged equipment
 * Prob   - the probability of the outage

  Examples:
(Outage generator at bus 1 with Id 1, probability = 0.1)
    1,0,1,0,0,1 ,0,0.1
(Outage branch connecting buses 8-9 with Id 1, probability = 0.1)
    2,1,0,8,9,1 ,0,0.1
(Multiple outages
  generator at bus 1 with Id 1, probability = 0.1
  branch connecting buses 8-9 with Id 1, probability = 0.1)
    3,0,1,0,0,1 ,0,0.1
    3,1,0,8,9,1 ,0,0.1
 */
PETSC_EXTERN PetscErrorCode
SCOPFLOWSetContingencyData(SCOPFLOW, ContingencyFileInputFormat, const char[]);
PETSC_EXTERN PetscErrorCode SCOPFLOWSetPLoadData(SCOPFLOW, const char[]);
PETSC_EXTERN PetscErrorCode SCOPFLOWSetQLoadData(SCOPFLOW, const char[]);
PETSC_EXTERN PetscErrorCode SCOPFLOWSetWindGenProfile(SCOPFLOW, const char[]);
PETSC_EXTERN PetscErrorCode SCOPFLOWSetTimeStep(SCOPFLOW, PetscReal);
PETSC_EXTERN PetscErrorCode SCOPFLOWSetDuration(SCOPFLOW, PetscReal);
PETSC_EXTERN PetscErrorCode SCOPFLOWSetUp(SCOPFLOW);
PETSC_EXTERN PetscErrorCode SCOPFLOWCreateGlobalVector(SCOPFLOW, Vec *);
PETSC_EXTERN PetscErrorCode SCOPFLOWCreateMatrix(SCOPFLOW, Mat *);
PETSC_EXTERN PetscErrorCode SCOPFLOWSolve(SCOPFLOW);
PETSC_EXTERN PetscErrorCode SCOPFLOWSetNumContingencies(SCOPFLOW, PetscInt);
PETSC_EXTERN PetscErrorCode SCOPFLOWGetTotalObjective(SCOPFLOW, PetscReal *);
PETSC_EXTERN PetscErrorCode SCOPFLOWGetBaseObjective(SCOPFLOW, PetscReal *);
PETSC_EXTERN PetscErrorCode SCOPFLOWGetSolution(SCOPFLOW, PetscInt, Vec *);
PETSC_EXTERN PetscErrorCode SCOPFLOWPrintSolution(SCOPFLOW, PetscInt);
PETSC_EXTERN PetscErrorCode SCOPFLOWSaveSolution(SCOPFLOW, PetscInt,
                                                 OutputFormat, const char[]);
PETSC_EXTERN PetscErrorCode SCOPFLOWSaveSolutionAll(SCOPFLOW, OutputFormat,
                                                    const char[]);
PETSC_EXTERN PetscErrorCode SCOPFLOWGetConvergenceStatus(SCOPFLOW, PetscBool *);
PETSC_EXTERN PetscErrorCode SCOPFLOWSetMode(SCOPFLOW, PetscInt);
PETSC_EXTERN PetscErrorCode SCOPFLOWGetNumIterations(SCOPFLOW, PetscInt *);
PETSC_EXTERN PetscErrorCode SCOPFLOWSetTolerance(SCOPFLOW, PetscReal);
PETSC_EXTERN PetscErrorCode SCOPFLOWGetTolerance(SCOPFLOW, PetscReal *);

PETSC_EXTERN PetscErrorCode SCOPFLOWSetTimeStepandDuration(SCOPFLOW, PetscReal,
                                                           PetscReal);
PETSC_EXTERN PetscErrorCode SCOPFLOWSetLoadProfiles(SCOPFLOW, const char[],
                                                    const char[]);
PETSC_EXTERN PetscErrorCode SCOPFLOWSetSubproblemModel(SCOPFLOW, const char[]);
PETSC_EXTERN PetscErrorCode SCOPFLOWSetSubproblemSolver(SCOPFLOW, const char[]);
PETSC_EXTERN PetscErrorCode SCOPFLOWSetComputeMode(SCOPFLOW, const char[]);
PETSC_EXTERN PetscErrorCode SCOPFLOWSetVerbosityLevel(SCOPFLOW, PetscInt);
PETSC_EXTERN PetscErrorCode SCOPFLOWEnablePowerImbalanceVariables(SCOPFLOW,
                                                                  PetscBool);
PETSC_EXTERN PetscErrorCode SCOPFLOWIgnoreLineflowConstraints(SCOPFLOW,
                                                              PetscBool);

#endif
