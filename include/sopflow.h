/*
 * Public header file for stochastic optimal power flow application.
 *
 */

#ifndef SOPFLOW_H
#define SOPFLOW_H

#include <opflow.h>
#include <scopflow.h>
#include <ps.h>

/* Models */
#define SOPFLOWMODEL_GENRAMP "GENRAMP"
/* Model for multi-contingency SOPFLOW */
#define SOPFLOWMODEL_GENRAMPC "GENRAMPC"

/* Solvers */
#define SOPFLOWSOLVER_IPOPT "IPOPT"
/* Embarassingly parallel solver - solves each OPFLOW independently */
#define SOPFLOWSOLVER_EMPAR "EMPAR"
/* Primal decomposition-basedd HIOp solver */
#define SOPFLOWSOLVER_HIOP "HIOP"

/* Initialization and Parameters*/
#define SOPFLOW_INITIALIZATION "ACPF"
#define SOPFLOW_GENBUSVOLTAGE "VARIABLE_WITHIN_BOUNDS"

/* Type of uncertainty */
typedef enum { NONE, WIND, LOAD } ScenarioUncertaintyType;

typedef struct _p_SOPFLOW *SOPFLOW;

namespace SOPFLOWOptions {

const auto sopflow_model = ExaGOStringOption(
    "-sopflow_model", "SOPFLOW model type", "GENRAMP", {"GENRAMPC"});

const auto sopflow_solver =
    ExaGOStringOption("-sopflow_solver", "SOPFLOW solver type",
                      SCOPFLOWOptions::solver.default_value,
                      SCOPFLOWOptions::solver.possible_values);

/* Retain default solver and model values from OPFLOW */

const auto opflow_model = OPFLOWOptions::model;

const auto subproblem_model = ExaGOStringOption(
    "-sopflow_subproblem_model", "SOPFLOW subproblem model type",
    OPFLOWOptions::model.default_value, OPFLOWOptions::model.possible_values);

const auto subproblem_solver = ExaGOStringOption(
    "-sopflow_subproblem_solver", "SOPFLOW subproblem solver type",
    OPFLOWOptions::solver.default_value, OPFLOWOptions::solver.possible_values);

const auto iscoupling = ExaGOBoolOption(
    "-sopflow_iscoupling",
    "Include coupling between first stage and second stage", PETSC_TRUE);

const auto Ns = ExaGOIntOption("-sopflow_Ns", "Number of scenarios", 1);
const auto Nc = ExaGOIntOption(
    "-sopflow_Nc", "Number of contingencies for multi-contingency scenario", 0);

const auto mode = ExaGOIntOption(
    "-sopflow_mode", "Operation mode:Preventive (0) or Corrective (1)", 1);
const auto enable_multicontingency =
    ExaGOBoolOption("-sopflow_enable_multicontingency",
                    "Multi-contingency SOPFLOW?", PETSC_FALSE);
const auto flatten_contingencies =
    ExaGOBoolOption("-sopflow_flatten_contingencies",
                    "Flatten contingencies for SOPFLOW?", PETSC_FALSE);

const auto tolerance =
    ExaGORealOption("-sopflow_tolerance", "Optimization tolerance", 1e-6);

const auto ctgcfile = ExaGOStringOption("-ctgcfile", "Contingency file",
                                        "/path/to/contingency_file", {});

const auto windgen = ExaGOStringOption("-windgen", "Wind generation file",
                                       "/path/to/windgen_file", {});

} // namespace SOPFLOWOptions

PETSC_EXTERN PetscErrorCode SOPFLOWEnableMultiContingency(SOPFLOW, PetscBool);
PETSC_EXTERN PetscErrorCode SOPFLOWSetModel(SOPFLOW, const char[]);
PETSC_EXTERN PetscErrorCode
SOPFLOWModelRegister(SOPFLOW, const char[], PetscErrorCode (*create)(SOPFLOW));
PETSC_EXTERN PetscErrorCode SOPFLOWSetSolver(SOPFLOW, const char[]);
PETSC_EXTERN PetscErrorCode SOPFLOWCreate(MPI_Comm, SOPFLOW *);
PETSC_EXTERN PetscErrorCode
    SOPFLOWSetInitializationType(SOPFLOW, OPFLOWInitializationType);
PETSC_EXTERN PetscErrorCode
    SOPFLOWSetGenBusVoltageType(SOPFLOW, OPFLOWGenBusVoltageType);
PETSC_EXTERN PetscErrorCode SOPFLOWDestroy(SOPFLOW *);
PETSC_EXTERN PetscErrorCode SOPFLOWSetNetworkData(SOPFLOW, const char[]);
PETSC_EXTERN PetscErrorCode SOPFLOWSetWindGenProfile(SOPFLOW, const char[]);

PETSC_EXTERN PetscErrorCode SOPFLOWSetScenarioData(SOPFLOW,
                                                   ScenarioFileInputFormat,
                                                   ScenarioUncertaintyType,
                                                   const char[]);
PETSC_EXTERN PetscErrorCode
SOPFLOWSetContingencyData(SOPFLOW, ContingencyFileInputFormat, const char[]);
PETSC_EXTERN PetscErrorCode SOPFLOWSetUp(SOPFLOW);
PETSC_EXTERN PetscErrorCode SOPFLOWCreateGlobalVector(SOPFLOW, Vec *);
PETSC_EXTERN PetscErrorCode SOPFLOWCreateMatrix(SOPFLOW, Mat *);
PETSC_EXTERN PetscErrorCode SOPFLOWSolve(SOPFLOW);
PETSC_EXTERN PetscErrorCode SOPFLOWSetNumScenarios(SOPFLOW, PetscInt);
PETSC_EXTERN PetscErrorCode SOPFLOWGetNumScenarios(SOPFLOW,
                                                   ScenarioFileInputFormat,
                                                   const char *, PetscInt *);
PETSC_EXTERN PetscErrorCode SOPFLOWFlattenContingencies(SOPFLOW, PetscBool);
PETSC_EXTERN PetscErrorCode SOPFLOWSetNumContingencies(SOPFLOW, PetscInt);
PETSC_EXTERN PetscErrorCode SOPFLOWGetTotalObjective(SOPFLOW, PetscReal *);
PETSC_EXTERN PetscErrorCode SOPFLOWGetBaseObjective(SOPFLOW, PetscReal *);
PETSC_EXTERN PetscErrorCode SOPFLOWGetSolution(SOPFLOW, PetscInt, Vec *);
PETSC_EXTERN PetscErrorCode SOPFLOWPrintSolution(SOPFLOW, PetscInt);
PETSC_EXTERN PetscErrorCode SOPFLOWSaveSolution(SOPFLOW, PetscInt, OutputFormat,
                                                const char[]);
PETSC_EXTERN PetscErrorCode SOPFLOWSaveSolutionAll(SOPFLOW, OutputFormat,
                                                   const char[]);
PETSC_EXTERN PetscErrorCode SOPFLOWGetConvergenceStatus(SOPFLOW, PetscBool *);
PETSC_EXTERN PetscErrorCode SOPFLOWGetMode(SOPFLOW, PetscInt *);
PETSC_EXTERN PetscErrorCode SOPFLOWGetNumIterations(SOPFLOW, PetscInt *);
PETSC_EXTERN PetscErrorCode SOPFLOWSetTolerance(SOPFLOW, PetscReal);
PETSC_EXTERN PetscErrorCode SOPFLOWGetTolerance(SOPFLOW, PetscReal *);

PETSC_EXTERN PetscErrorCode SOPFLOWSetTimeStepandDuration(SOPFLOW, PetscReal,
                                                          PetscReal);
PETSC_EXTERN PetscErrorCode SOPFLOWSetLoadProfiles(SOPFLOW, const char[],
                                                   const char[]);
PETSC_EXTERN PetscErrorCode SOPFLOWSetSubproblemModel(SOPFLOW, const char[]);
PETSC_EXTERN PetscErrorCode SOPFLOWSetSubproblemSolver(SOPFLOW, const char[]);
PETSC_EXTERN PetscErrorCode SOPFLOWSetIgnoreLineflowConstraints(SOPFLOW,
                                                                PetscBool);
PETSC_EXTERN PetscErrorCode SOPFLOWSetSubproblemComputeMode(SOPFLOW,
                                                            const char[]);
PETSC_EXTERN PetscErrorCode SOPFLOWSetSubproblemVerbosityLevel(SOPFLOW, int);

#endif
