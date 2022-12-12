#include "exago_python_sopflow.hpp"

// -------------------------------------------------------------
// class SOPFLOW_wrapper
//
// Wrap to make sure allocation and destruction is handled correctly
// -------------------------------------------------------------
class SOPFLOW_wrapper {
public:
  SOPFLOW sopf;
  SOPFLOW_wrapper(void) {
    PetscErrorCode ierr;
    MPI_Comm communicator;
    ierr = ExaGOGetSelfCommunicator(&communicator);
    ExaGOCheckError(ierr);
    ierr = SOPFLOWCreate(communicator, &sopf);
    ExaGOCheckError(ierr);
  }
  ~SOPFLOW_wrapper(void) {
    PetscErrorCode ierr;
    ierr = SOPFLOWDestroy(&sopf);
    ExaGOCheckError(ierr);
  }
};

// -------------------------------------------------------------
// init_exago_sopflow
// -------------------------------------------------------------
void init_exago_sopflow(pybind11::module &m) {

  pybind11::enum_<ScenarioFileInputFormat>(m, "ScenarioFileInputFormat")
      .value("NATIVE_SINGLEPERIOD", SOPFLOW_NATIVE_SINGLEPERIOD)
      .value("NATIVE_MULTIPERIOD", SOPFLOW_NATIVE_MULTIPERIOD)
      .export_values();

  pybind11::enum_<ScenarioUncertaintyType>(m, "ScenarioUncertaintyType")
      .value("NONE", NONE)
      .value("WIND", WIND)
      .value("LOAD", LOAD)
      .export_values();

  pybind11::class_<SOPFLOW_wrapper>(m, "SOPFLOW")
      .def(pybind11::init())

      /* Setters */

      .def("set_model",
           [](SOPFLOW_wrapper &w, std::string model) {
             PetscErrorCode ierr;
             ierr = SOPFLOWSetModel(w.sopf, model.c_str());
             ExaGOCheckError(ierr);
           })
      .def("set_network_data",
           [](SOPFLOW_wrapper &w, std::string filename) {
             PetscErrorCode ierr;
             ierr = SOPFLOWSetNetworkData(w.sopf, filename.c_str());
             ExaGOCheckError(ierr);
           })
      .def("set_contingency_data",
           [](SOPFLOW_wrapper &w, std::string filename,
              ContingencyFileInputFormat fmt = NATIVE) {
             PetscErrorCode ierr;
             ierr = SOPFLOWSetContingencyData(w.sopf, fmt, filename.c_str());
             ExaGOCheckError(ierr);
           })
      .def("set_num_contingencies",
           [](SOPFLOW_wrapper &w, int n) {
             PetscErrorCode ierr;
             ierr = SOPFLOWSetNumContingencies(w.sopf, n);
             ExaGOCheckError(ierr);
           })
      .def("set_scenario_data",
           [](SOPFLOW_wrapper &w, std::string filename,
              ScenarioFileInputFormat sfmt, ScenarioUncertaintyType ut) {
             PetscErrorCode ierr;
             ierr = SOPFLOWSetScenarioData(w.sopf, sfmt, ut, filename.c_str());
             ExaGOCheckError(ierr);
           })
      .def("set_num_scenarios",
           [](SOPFLOW_wrapper &w, int n) {
             PetscErrorCode ierr;
             ierr = SOPFLOWSetNumScenarios(w.sopf, n);
             ExaGOCheckError(ierr);
           })
      .def("set_wind_gen_profile",
           [](SOPFLOW_wrapper &w, std::string filename) {
             PetscErrorCode ierr;
             ierr = SOPFLOWSetWindGenProfile(w.sopf, filename.c_str());
             ExaGOCheckError(ierr);
           })
      .def("set_time_step_and_duration",
           [](SOPFLOW_wrapper &w, double dt, double tmax) {
             PetscErrorCode ierr;
             ierr = SOPFLOWSetTimeStepandDuration(w.sopf, dt, tmax);
             ExaGOCheckError(ierr);
           })
      .def("set_tolerance",
           [](SOPFLOW_wrapper &w, double tol) {
             PetscErrorCode ierr;
             ierr = SOPFLOWSetTolerance(w.sopf, tol);
             ExaGOCheckError(ierr);
           })
      .def("set_subproblem_verbosity_level",
           [](SOPFLOW_wrapper &w, int level) {
             PetscErrorCode ierr;
             ierr = SOPFLOWSetSubproblemVerbosityLevel(w.sopf, level);
             ExaGOCheckError(ierr);
           })
      .def("set_subproblem_compute_mode",
           [](SOPFLOW_wrapper &w, std::string mode) {
             PetscErrorCode ierr;
             ierr = SOPFLOWSetSubproblemComputeMode(w.sopf, mode.c_str());
             ExaGOCheckError(ierr);
           })
      .def("set_subproblem_model",
           [](SOPFLOW_wrapper &w, std::string model) {
             PetscErrorCode ierr;
             ierr = SOPFLOWSetSubproblemModel(w.sopf, model.c_str());
             ExaGOCheckError(ierr);
           })
      .def("set_subproblem_solver",
           [](SOPFLOW_wrapper &w, std::string solver) {
             PetscErrorCode ierr;
             ierr = SOPFLOWSetSubproblemSolver(w.sopf, solver.c_str());
             ExaGOCheckError(ierr);
           })
      .def("set_solver",
           [](SOPFLOW_wrapper &w, std::string solver) {
             PetscErrorCode ierr;
             ierr = SOPFLOWSetSolver(w.sopf, solver.c_str());
             ExaGOCheckError(ierr);
           })

      .def("set_initialization_type",
           [](SOPFLOW_wrapper &w, std::string init) {
             PetscErrorCode ierr;
             OPFLOWInitializationType type;
             bool found = false;
#define OPFLOW_INIT_DEF(x)                                                     \
  if (init == #x) {                                                            \
    found = true;                                                              \
    type = x;                                                                  \
  }
#include "private/opflow_enum.def"
#undef OPFLOW_INIT_DEF
             if (!found) {
#ifdef EXAGO_ENABLE_LOGGING
               throw ExaGOError(
                   fmt::format(
                       "Invalid value '{}' for OPFLOWInitializationType enum",
                       init)
                       .c_str());
#else
             char sbuf[256];
             sprintf(sbuf,"Invalid value %s for OPFLOWInitializationType enum",init.c_str());
             throw ExaGOError(sbuf);
#endif
             }
             ierr = SOPFLOWSetInitializationType(w.sopf, type);
             ExaGOCheckError(ierr);
           })
      .def("set_initialization_type",
           [](SOPFLOW_wrapper &w, int init) {
             PetscErrorCode ierr;
             OPFLOWInitializationType type;
             bool found = false;
#define OPFLOW_INIT_DEF(x)                                                     \
  if (init == x) {                                                             \
    found = true;                                                              \
    type = x;                                                                  \
  }
#include "private/opflow_enum.def"
#undef OPFLOW_INIT_DEF
             if (!found) {
#ifdef EXAGO_ENABLE_LOGGING
               throw ExaGOError(
                   fmt::format(
                       "Invalid value '{}' for SOPFLOWInitializationType enum",
                       init)
                       .c_str());
#else
             char sbuf[256];
             sprintf(sbuf,"Invalid value %d for SOPFLOWInitializationType enum",init);
             throw ExaGOError(sbuf);
#endif
             }
             ierr = SOPFLOWSetInitializationType(w.sopf, type);
             ExaGOCheckError(ierr);
           })
      .def("set_gen_bus_voltage_type",
           [](SOPFLOW_wrapper &w, int volt) {
             PetscErrorCode ierr;
             OPFLOWGenBusVoltageType type;
             bool found = false;
#define OPFLOW_GBV_DEF(x)                                                      \
  if (volt == x) {                                                             \
    found = true;                                                              \
    type = x;                                                                  \
  }
#include "private/opflow_enum.def"
#undef OPFLOW_GBV_DEF
             if (!found) {
#ifdef EXAGO_ENABLE_LOGGING
               throw ExaGOError(
                   fmt::format(
                       "Invalid value '{}' for SOPFLOWGenBusVoltageType enum",
                       volt)
                       .c_str());
#else
             char sbuf[256];
             sprintf(sbuf,"Invalid value %d for SOPFLOWGenBusVoltageType enum",volt);
             throw ExaGOError(sbuf);
#endif
             }
             ierr = SOPFLOWSetGenBusVoltageType(w.sopf, type);
             ExaGOCheckError(ierr);
           })

      .def("set_gen_bus_voltage_type",
           [](SOPFLOW_wrapper &w, std::string volt) {
             PetscErrorCode ierr;
             OPFLOWGenBusVoltageType type;
             bool found = false;
#define OPFLOW_GBV_DEF(x)                                                      \
  if (volt == #x) {                                                            \
    found = true;                                                              \
    type = x;                                                                  \
  }
#include "private/opflow_enum.def"
#undef OPFLOW_GBV_DEF
             if (!found) {
#ifdef EXAGO_ENABLE_LOGGING
               throw ExaGOError(
                   fmt::format(
                       "Invalid value '{}' for SOPFLOWGenBusVoltageType enum",
                       volt)
                       .c_str());
#else
             char sbuf[256];
             sprintf(sbuf,"Invalid value %s for SOPFLOWGenBusVoltageType enum",volt.c_str());
             throw ExaGOError(sbuf);
#endif
             }
             ierr = SOPFLOWSetGenBusVoltageType(w.sopf, type);
             ExaGOCheckError(ierr);
           })
      .def("enable_multi_contingency",
           [](SOPFLOW_wrapper &w, bool flag) {
             PetscErrorCode ierr;
             ierr = SOPFLOWEnableMultiContingency(w.sopf, (PetscBool)flag);
             ExaGOCheckError(ierr);
           })

      /* Getters */

      .def("get_num_scenarios",
           [](SOPFLOW_wrapper &w, const std::string &filename,
              ScenarioFileInputFormat fmt) -> int {
             PetscErrorCode ierr;
             PetscInt n;
             ierr = SOPFLOWGetNumScenarios(w.sopf, fmt, filename.c_str(), &n);
             ExaGOCheckError(ierr);
             return n;
           })
      .def("get_num_iterations",
           [](SOPFLOW_wrapper &w) -> int {
             PetscErrorCode ierr;
             PetscInt n;
             ierr = SOPFLOWGetNumIterations(w.sopf, &n);
             ExaGOCheckError(ierr);
             return n;
           })
      .def("get_convergence_status",
           [](SOPFLOW_wrapper &w) -> bool {
             PetscErrorCode ierr;
             PetscBool flag;
             ierr = SOPFLOWGetConvergenceStatus(w.sopf, &flag);
             ExaGOCheckError(ierr);
             return flag;
           })
      .def("get_total_objective",
           [](SOPFLOW_wrapper &w) -> double {
             PetscErrorCode ierr;
             double obj;
             ierr = SOPFLOWGetTotalObjective(w.sopf, &obj);
             ExaGOCheckError(ierr);
             return obj;
           })
      .def("get_base_objective",
           [](SOPFLOW_wrapper &w) -> double {
             PetscErrorCode ierr;
             double obj;
             ierr = SOPFLOWGetBaseObjective(w.sopf, &obj);
             ExaGOCheckError(ierr);
             return obj;
           })
      .def("get_converged_status",
           [](SOPFLOW_wrapper &w) -> bool {
             PetscErrorCode ierr;
             PetscBool flag;
             ierr = SOPFLOWGetConvergenceStatus(w.sopf, &flag);
             return (bool)flag;
           })
      .def("get_tolerance",
           [](SOPFLOW_wrapper &w) -> double {
             PetscErrorCode ierr;
             double tol;
             ierr = SOPFLOWGetTolerance(w.sopf, &tol);
             ExaGOCheckError(ierr);
             return tol;
           })

      /* Solution */

      .def("setup",
           [](SOPFLOW_wrapper &w) {
             PetscErrorCode ierr;
             ierr = SOPFLOWSetUp(w.sopf);
             ExaGOCheckError(ierr);
           })
      .def("solve",
           [](SOPFLOW_wrapper &w) {
             PetscErrorCode ierr;
             ierr = SOPFLOWSolve(w.sopf);
             ExaGOCheckError(ierr);
           })

      ;
}
