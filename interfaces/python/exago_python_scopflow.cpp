#include "exago_python_scopflow.hpp"

// -------------------------------------------------------------
// class SCOPFLOW_wrapper
//
// Wrap to make sure allocation and destruction is handled correctly
// -------------------------------------------------------------
class SCOPFLOW_wrapper {
public:
  SCOPFLOW scopf;
  SCOPFLOW_wrapper(void) {
    PetscErrorCode ierr;
    MPI_Comm communicator;
    ierr = ExaGOGetSelfCommunicator(&communicator);
    ExaGOCheckError(ierr);
    ierr = SCOPFLOWCreate(communicator, &scopf);
    ExaGOCheckError(ierr);
  }
  ~SCOPFLOW_wrapper(void) {
    PetscErrorCode ierr;
    ierr = SCOPFLOWDestroy(&scopf);
    ExaGOCheckError(ierr);
  }
};

// -------------------------------------------------------------
// init_exago_scopflow
// -------------------------------------------------------------
void init_exago_scopflow(pybind11::module &m) {

  pybind11::enum_<ContingencyFileInputFormat>(m, "ContingencyFileInputFormat")
      .value("NATIVE", NATIVE)
      .value("PSSE", PSSE)
      .export_values();

  pybind11::class_<SCOPFLOW_wrapper>(m, "SCOPFLOW")
      .def(pybind11::init())
      /* Setters */
      .def("set_model",
           [](SCOPFLOW_wrapper &w, std::string model) {
             PetscErrorCode ierr;
             ierr = SCOPFLOWSetModel(w.scopf, model.c_str());
             ExaGOCheckError(ierr);
           })
      .def("set_network_data",
           [](SCOPFLOW_wrapper &w, std::string filename) {
             PetscErrorCode ierr;
             ierr = SCOPFLOWSetNetworkData(w.scopf, filename.c_str());
             ExaGOCheckError(ierr);
           })
      .def("set_num_contingencies",
           [](SCOPFLOW_wrapper &w, int n) {
             PetscErrorCode ierr;
             ierr = SCOPFLOWSetNumContingencies(w.scopf, n);
             ExaGOCheckError(ierr);
           })
      .def("set_contingency_data",
           [](SCOPFLOW_wrapper &w, std::string filename,
              ContingencyFileInputFormat fmt) {
             PetscErrorCode ierr;
             ierr = SCOPFLOWSetContingencyData(w.scopf, fmt, filename.c_str());
             ExaGOCheckError(ierr);
           })
      .def("set_pload_data",
           [](SCOPFLOW_wrapper &w, std::string filename) {
             PetscErrorCode ierr;
             ierr = SCOPFLOWSetPLoadData(w.scopf, filename.c_str());
             ExaGOCheckError(ierr);
           })
      .def("set_qload_data",
           [](SCOPFLOW_wrapper &w, std::string filename) {
             PetscErrorCode ierr;
             ierr = SCOPFLOWSetQLoadData(w.scopf, filename.c_str());
             ExaGOCheckError(ierr);
           })
      .def("set_load_profiles",
           [](SCOPFLOW_wrapper &w, const std::string &pfile,
              const std::string &qfile) {
             PetscErrorCode ierr;
             ierr =
                 SCOPFLOWSetLoadProfiles(w.scopf, pfile.c_str(), qfile.c_str());
             ExaGOCheckError(ierr);
           })
      .def("set_wind_gen_profile",
           [](SCOPFLOW_wrapper &w, std::string filename) {
             PetscErrorCode ierr;
             ierr = SCOPFLOWSetWindGenProfile(w.scopf, filename.c_str());
             ExaGOCheckError(ierr);
           })
      .def("set_time_step",
           [](SCOPFLOW_wrapper &w, double dt) {
             PetscErrorCode ierr;
             ierr = SCOPFLOWSetTimeStep(w.scopf, dt);
             ExaGOCheckError(ierr);
           })
      .def("set_duration",
           [](SCOPFLOW_wrapper &w, double d) {
             PetscErrorCode ierr;
             ierr = SCOPFLOWSetDuration(w.scopf, d);
             ExaGOCheckError(ierr);
           })
      .def("set_time_step_and_duration",
           [](SCOPFLOW_wrapper &w, double dt, double d) {
             PetscErrorCode ierr;
             ierr = SCOPFLOWSetTimeStepandDuration(w.scopf, dt, d);
             ExaGOCheckError(ierr);
           })
      .def("set_tolerance",
           [](SCOPFLOW_wrapper &w, double tol) {
             PetscErrorCode ierr;
             ierr = SCOPFLOWSetTolerance(w.scopf, tol);
             ExaGOCheckError(ierr);
           })
      .def("set_verbosity_level",
           [](SCOPFLOW_wrapper &w, int level) {
             PetscErrorCode ierr;
             ierr = SCOPFLOWSetSubproblemVerbosityLevel(w.scopf, level);
             ExaGOCheckError(ierr);
           })
      .def("set_compute_mode",
           [](SCOPFLOW_wrapper &w, std::string mode) {
             PetscErrorCode ierr;
             ierr = SCOPFLOWSetSubproblemComputeMode(w.scopf, mode.c_str());
             ExaGOCheckError(ierr);
           })
      .def("set_mode",
           [](SCOPFLOW_wrapper &w, int mode) {
             PetscErrorCode ierr;
             ierr = SCOPFLOWSetMode(w.scopf, mode);
             ExaGOCheckError(ierr);
           })
      .def("set_solver",
           [](SCOPFLOW_wrapper &w, std::string solver) {
             PetscErrorCode ierr;
             ierr = SCOPFLOWSetSolver(w.scopf, solver.c_str());
             ExaGOCheckError(ierr);
           })
      .def("set_subproblem_model",
           [](SCOPFLOW_wrapper &w, std::string model) {
             PetscErrorCode ierr;
             ierr = SCOPFLOWSetSubproblemModel(w.scopf, model.c_str());
             ExaGOCheckError(ierr);
           })
      .def("set_subproblem_solver",
           [](SCOPFLOW_wrapper &w, std::string solver) {
             PetscErrorCode ierr;
             ierr = SCOPFLOWSetSubproblemSolver(w.scopf, solver.c_str());
             ExaGOCheckError(ierr);
           })
      .def("set_initialization_type",
           [](SCOPFLOW_wrapper &w, std::string init) {
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
             ierr = SCOPFLOWSetInitilizationType(w.scopf, type);
             ExaGOCheckError(ierr);
           })
      .def("set_initialization_type",
           [](SCOPFLOW_wrapper &w, int init) {
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
                       "Invalid value '{}' for SCOPFLOWInitializationType enum",
                       init)
                       .c_str());
#else
             char sbuf[256];
             sprintf(sbuf,"Invalid value %d for OPFLOWInitializationType enum",init);
             throw ExaGOError(sbuf);
#endif
             }
             ierr = SCOPFLOWSetInitilizationType(w.scopf, type);
             ExaGOCheckError(ierr);
           })

      .def("set_gen_bus_voltage_type",
           [](SCOPFLOW_wrapper &w, int volt) {
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
                       "Invalid value '{}' for SCOPFLOWGenBusVoltageType enum",
                       volt)
                       .c_str());
#else
             char sbuf[256];
             sprintf(sbuf,"Invalid value %d for OPFLOWGenBusVoltageType enum",volt);
             throw ExaGOError(sbuf);
#endif
             }
             ierr = SCOPFLOWSetGenBusVoltageType(w.scopf, type);
             ExaGOCheckError(ierr);
           })

      .def("set_gen_bus_voltage_type",
           [](SCOPFLOW_wrapper &w, std::string volt) {
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
                       "Invalid value '{}' for SCOPFLOWGenBusVoltageType enum",
                       volt)
                       .c_str());
#else
             char sbuf[256];
             sprintf(sbuf,"Invalid value %s for OPFLOWGenBusVoltageType enum",volt.c_str());
             throw ExaGOError(sbuf);
#endif
             }
             ierr = SCOPFLOWSetGenBusVoltageType(w.scopf, type);
             ExaGOCheckError(ierr);
           })

      .def("enable_multi_period",
           [](SCOPFLOW_wrapper &w, bool flag) {
             PetscErrorCode ierr;
             ierr = SCOPFLOWEnableMultiPeriod(w.scopf, (PetscBool)flag);
             ExaGOCheckError(ierr);
           })
      .def("enable_power_imbalance_variables",
           [](SCOPFLOW_wrapper &w, bool flag) {
             PetscErrorCode ierr;
             ierr = SCOPFLOWEnablePowerImbalanceVariables(w.scopf,
                                                          (PetscBool)flag);
             ExaGOCheckError(ierr);
           })
      .def("ignore_lineflow_constraints",
           [](SCOPFLOW_wrapper &w, bool flag) {
             PetscErrorCode ierr;
             ierr = SCOPFLOWIgnoreLineflowConstraints(w.scopf, (PetscBool)flag);
             ExaGOCheckError(ierr);
           })

      /* Getters */
      .def("get_tolerance",
           [](SCOPFLOW_wrapper &w) -> double {
             PetscErrorCode ierr;
             PetscReal tol;
             ierr = SCOPFLOWGetTolerance(w.scopf, &tol);
             ExaGOCheckError(ierr);
             return tol;
           })
      .def("get_num_iterations",
           [](SCOPFLOW_wrapper &w) -> int {
             PetscErrorCode ierr;
             PetscInt n;
             ierr = SCOPFLOWGetNumIterations(w.scopf, &n);
             ExaGOCheckError(ierr);
             return n;
           })
      .def("get_convergence_status",
           [](SCOPFLOW_wrapper &w) -> bool {
             PetscErrorCode ierr;
             PetscBool flag;
             ierr = SCOPFLOWGetConvergenceStatus(w.scopf, &flag);
             ExaGOCheckError(ierr);
             return flag;
           })
      .def("get_total_objective",
           [](SCOPFLOW_wrapper &w) -> double {
             PetscErrorCode ierr;
             PetscReal obj;
             ierr = SCOPFLOWGetTotalObjective(w.scopf, &obj);
             ExaGOCheckError(ierr);
             return obj;
           })
      .def("get_base_objective",
           [](SCOPFLOW_wrapper &w) -> double {
             PetscErrorCode ierr;
             PetscReal obj;
             ierr = SCOPFLOWGetBaseObjective(w.scopf, &obj);
             ExaGOCheckError(ierr);
             return obj;
           })

      /* Solution */
      .def("set_up",
           [](SCOPFLOW_wrapper &w) {
             PetscErrorCode ierr;
             ierr = SCOPFLOWSetUp(w.scopf);
             ExaGOCheckError(ierr);
           })
      .def("solve",
           [](SCOPFLOW_wrapper &w) {
             PetscErrorCode ierr;
             ierr = SCOPFLOWSolve(w.scopf);
             ExaGOCheckError(ierr);
           })

      .def("print_solution",
           [](SCOPFLOW_wrapper &w, int n) {
             PetscErrorCode ierr;
             ierr = SCOPFLOWPrintSolution(w.scopf, n);
             ExaGOCheckError(ierr);
           })

      .def("save_solution",
           [](SCOPFLOW_wrapper &w, int n, OutputFormat fmt,
              std::string outfile) {
             PetscErrorCode ierr;
             ierr = SCOPFLOWSaveSolution(w.scopf, n, fmt, outfile.c_str());
             ExaGOCheckError(ierr);
           })

      .def("save_solution_all",
           [](SCOPFLOW_wrapper &w, OutputFormat fmt, std::string outdir) {
             PetscErrorCode ierr;
             ierr = SCOPFLOWSaveSolutionAll(w.scopf, fmt, outdir.c_str());
             ExaGOCheckError(ierr);
           });
}
